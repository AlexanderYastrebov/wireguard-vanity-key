// Package main searches for vanity X25519 key based on algorithm
// described here https://github.com/warner/wireguard-vanity-address/pull/15
package main

import (
	"bytes"
	"context"
	"crypto/rand"
	"encoding/base64"
	"encoding/binary"
	"flag"
	"fmt"
	"io"
	"math/big"
	"os"
	"os/signal"
	"runtime"
	"sync"
	"sync/atomic"
	"syscall"
	"time"

	"filippo.io/edwards25519"
	"filippo.io/edwards25519/field"
)

var (
	// l = 2^252 + 27742317777372353535851937790883648493 == 2^252 + smallScalar
	smallScalar = scalarFromBytes(decimalToBytes("27742317777372353535851937790883648493")...)
	// Ed25519 group's cofactor
	scalarOffset = scalarFromBytes(8)
	pointOffset  = new(edwards25519.Point).ScalarBaseMult(scalarOffset)
)

func main() {
	start := time.Now()

	prefix := flag.String("prefix", "AY/", "prefix of base64-encoded public key")
	timeout := flag.Duration("timeout", 0, "stop after specified timeout")
	flag.Parse()

	s0, p0 := newPair()

	ctx, cancel := signal.NotifyContext(context.Background(), syscall.SIGINT, syscall.SIGTERM)
	defer cancel()

	if *timeout != 0 {
		ctx, cancel = context.WithTimeout(ctx, *timeout)
		defer cancel()
	}

	test := func(p []byte) bool {
		return hasBase64Prefix(p, []byte(*prefix))
	}

	p, n, attempts := findPointParallel(ctx, runtime.NumCPU(), p0, test)

	private := "-"
	public := *prefix + "..."
	if p != nil {
		s := adjustScalar(s0, n)
		private = base64.StdEncoding.EncodeToString(scalarToKeyBytes(s))
		public = base64.StdEncoding.EncodeToString(p.BytesMontgomery())
	}

	duration := time.Since(start)

	fmt.Printf("%-44s %-44s %-10s %-10s %s\n", "private", "public", "attempts", "duration", "attempts/s")
	fmt.Printf("%-44s %-44s %-10d %-10s %d\n", private, public, attempts, duration.Round(time.Second), time.Duration(attempts)*(time.Second)/duration)

	if p == nil {
		os.Exit(1)
	}
}

func newPair() (*edwards25519.Scalar, *edwards25519.Point) {
	var key [32]byte
	if _, err := io.ReadFull(rand.Reader, key[:]); err != nil {
		panic(err)
	}

	s, err := edwards25519.NewScalar().SetBytesWithClamping(key[:])
	if err != nil {
		panic(err)
	}

	p := new(edwards25519.Point).ScalarBaseMult(s)

	return s, p
}

func findPointParallel(ctx context.Context, workers int, p0 *edwards25519.Point, test func([]byte) bool) (*edwards25519.Point, uint64, uint64) {
	type point struct {
		p *edwards25519.Point
		n uint64
	}

	result := make(chan point, workers)
	var attempts atomic.Uint64

	gctx, cancel := context.WithCancel(ctx)
	defer cancel()

	var wg sync.WaitGroup
	wg.Add(workers)
	for range workers {
		go func() {
			defer wg.Done()

			skip := randUint64()
			// p, n := findPoint(gctx, p0, skip, test)
			p, n := findBatchPoint(gctx, p0, skip, 1024, test)

			attempts.Add(n - skip)
			if p != nil {
				result <- point{p, n}
				cancel()
			}
		}()
	}
	wg.Wait()

	select {
	case r := <-result:
		return r.p, r.n, attempts.Load()
	case <-ctx.Done():
		return nil, 0, attempts.Load()
	}
}

func findPoint(ctx context.Context, p0 *edwards25519.Point, skip uint64, test func([]byte) bool) (*edwards25519.Point, uint64) {
	skipOffset := new(edwards25519.Point).ScalarMult(scalarFromUint64(skip), pointOffset)
	p := new(edwards25519.Point).Add(p0, skipOffset)
	n := skip

	var bm [32]byte
	bytesMontgomery(p, &bm)

	for ; !test(bm[:]); n++ {
		select {
		case <-ctx.Done():
			return nil, n
		default:
			p.Add(p, pointOffset)
			bytesMontgomery(p, &bm)
		}
	}
	return p, n
}

func findBatchPoint(ctx context.Context, p0 *edwards25519.Point, skip uint64, batchSize int, test func([]byte) bool) (*edwards25519.Point, uint64) {
	skipOffset := new(edwards25519.Point).ScalarMult(scalarFromUint64(skip), pointOffset)
	p := new(edwards25519.Point).Add(p0, skipOffset)

	n := skip

	pts := make([]edwards25519.Point, batchSize)
	u := make([]field.Element, batchSize)
	scratch := make([][]field.Element, 4)

	for i := range scratch {
		scratch[i] = make([]field.Element, batchSize)
	}

	var bm [32]byte
	for {
		select {
		case <-ctx.Done():
			return nil, n
		default:
		}

		for i := range pts {
			pts[i].Set(p)
			p.Add(p, pointOffset)
		}

		batchBytesMontgomery(pts, u, scratch)

		for i := range pts {
			copy(bm[:], u[i].Bytes()) // eliminate field.Element.Bytes() allocations
			if test(bm[:]) {
				return &pts[i], n + uint64(i)
			}
		}

		n += uint64(len(pts))
	}
}

func adjustScalar(s *edwards25519.Scalar, n uint64) *edwards25519.Scalar {
	return edwards25519.NewScalar().MultiplyAdd(scalarOffset, scalarFromUint64(n), s)
}

func scalarToKeyBytes(s *edwards25519.Scalar) []byte {
	// We can't use Scalars to add "l" and produce the aliases: any addition
	// we do on the Scalar will be reduced immediately. But we can add
	// "small", and then manually adjust the high-end byte, to produce an
	// array of bytes whose value is s+kl
	//
	// The aliases (with high probability) have distinct
	// high-order bits: 0b0001, 0b0010, etc. We want one of the four aliases
	// whose high-order bits are 0b01xx: these bits will survive the high-end
	// clamping unchanged. These are where k=[4..7].
	//
	// The three low-order bits will be some number N. Each alias adds l%8 == 5 to
	// this low end:
	//
	// $ echo '(2^252 + 27742317777372353535851937790883648493) % 8' | bc
	// 5
	//
	// So the first alias (k=1) will end in N+5, the second
	// (k=2) will end in N+2 (since (5+5)%8 == 2). Our k=4..7 yields
	// N+4,N+1,N+6,N+3. One of these values might be all zeros. That alias
	// will survive the low-end clamping unchanged.

	lowBits := s.Bytes()[0] & 0b111
	// Solve (lowBits + k*5) % 8 == 0 for k:
	// k := [8]byte{0, 0, 6, 0, 4, 7, 0, 5}[lowBits]
	k := [8]byte{0, 3, 6, 1, 4, 7, 2, 5}[lowBits]
	if k < 4 { // TODO: why k is mostly one of 4, 5, 6, 7 when scalarOffset is Ed25519 group's cofactor?
		panic("invalid scalar first byte (lowBits)")
	}

	aliasBytes := edwards25519.NewScalar().MultiplyAdd(smallScalar, scalarFromBytes(k), s).Bytes()
	aliasBytes[31] += (k << 4)

	return aliasBytes
}

func scalarFromBytes(x ...byte) *edwards25519.Scalar {
	var xb [64]byte
	copy(xb[:], x)

	xs, err := edwards25519.NewScalar().SetUniformBytes(xb[:])
	if err != nil {
		panic(err)
	}
	return xs
}

func scalarFromUint64(n uint64) *edwards25519.Scalar {
	var nb [8]byte
	binary.LittleEndian.PutUint64(nb[:], n)
	return scalarFromBytes(nb[:]...)
}

func decimalToBytes(d string) []byte {
	i, ok := new(big.Int).SetString(d, 10)
	if !ok {
		panic("invalid decimal string " + d)
	}
	b := i.Bytes()
	// convert to little endian
	for i, j := 0, len(b)-1; i < j; i, j = i+1, j-1 {
		b[i], b[j] = b[j], b[i]
	}
	return b
}

func hasBase64Prefix(p, prefix []byte) bool {
	var dst [44]byte
	base64.StdEncoding.Encode(dst[:], p)
	return bytes.HasPrefix(dst[:], prefix)
}

func randUint64() uint64 {
	var num uint64
	err := binary.Read(rand.Reader, binary.NativeEndian, &num)
	if err != nil {
		panic(err)
	}
	return num
}

// bytesMontgomery is a copy of [edwards25519.Point.BytesMontgomery]
// to eliminate allocations.
//
// bytesMontgomery uses:
//
//	1 addition
//	1 subtraction
//	1 invert = 254 squaring + 11 multiplications
//	1 multiplication
//
// i.e. ~254+11+1 = 266 multiplications
func bytesMontgomery(v *edwards25519.Point, buf *[32]byte) {
	// RFC 7748, Section 4.1 provides the bilinear map to calculate the
	// Montgomery u-coordinate
	//
	//              u = (1 + y) / (1 - y)
	//
	// where y = Y / Z and therefore
	//
	//              u = (Z + Y) / (Z - Y)

	var n, r, u field.Element

	_, Y, Z, _ := v.ExtendedCoordinates()
	n.Add(Z, Y)                // n = Z + Y
	r.Invert(r.Subtract(Z, Y)) // r = 1 / (Z - Y)
	u.Multiply(&n, &r)         // u = n * r

	copy(buf[:], u.Bytes())
}

// batchBytesMontgomery is equivalent to calling [edwards25519.Point.BytesMontgomery] for each point
// except that it uses [vectorDivision] and thus uses less point multiplications.
//
// All input slices must be of the same length n.
// Result bytes are encoded into u using scratch which should be at least 4 slices of length n.
//
// batchBytesMontgomery uses:
//
//	n additions
//	n subtractions
//	vectorDivision = 265+4*(n-1)+1 multiplications
//
// i.e. ~4*n multiplications for large n.
func batchBytesMontgomery(pts []edwards25519.Point, u []field.Element, scratch [][]field.Element) {
	x := scratch[0]
	y := scratch[1]

	// u = (Z + Y) / (Z - Y) = x / y
	for i, v := range pts {
		_, Y, Z, _ := v.ExtendedCoordinates()
		x[i].Add(Z, Y)      // x = Z + Y
		y[i].Subtract(Z, Y) // y = Z - Y
	}

	vectorDivision(x, y, u, scratch[2:]) // u = x / y
}

// vectorDivision calculates u = x / y using scratch.
//
// vectorDivision uses:
//
//	4*(n-1)+1 multiplications
//	1 invert = ~265 multiplications
//
// i.e. ~265+4*(n-1)+1 multiplications
//
// Simultaneous field divisions: an extension of Montgomery's trick
// David G. Harris
// https://eprint.iacr.org/2008/199.pdf
func vectorDivision(x, y, u []field.Element, scratch [][]field.Element) {
	n := len(x)
	r := scratch[0]
	s := scratch[1]

	r[0] = y[0]
	for i := 1; i < n; i++ {
		r[i].Multiply(&r[i-1], &y[i])
		s[i].Multiply(&r[i-1], &x[i])
	}

	I := new(field.Element).Invert(&r[n-1])

	t := I
	for i := n - 1; i > 0; i-- {
		u[i].Multiply(t, &s[i])
		t.Multiply(t, &y[i])
	}
	u[0].Multiply(t, &x[0])
}

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
	"strings"
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
	if len(os.Args) > 1 && os.Args[1] == "add" {
		cmdAdd(os.Args[2:])
		return
	}

	start := time.Now()
	config := struct {
		prefix  string
		timeout time.Duration
		public  string
		output  string
	}{}

	flag.StringVar(&config.prefix, "prefix", "AY/", "prefix of base64-encoded public key")
	flag.DurationVar(&config.timeout, "timeout", 0, "stop after specified timeout")
	flag.StringVar(&config.public, "public", "", "start from specified public key")
	flag.StringVar(&config.output, "output", "", "use \"offset\" to print offset only")
	flag.Parse()

	var s0 *edwards25519.Scalar
	var p0 *edwards25519.Point
	var err error

	if config.public != "" {
		p0, err = parsePublicKey(config.public)
		if err != nil {
			panic(err)
		}
	} else {
		s0, p0 = newPair()
	}

	ctx, cancel := signal.NotifyContext(context.Background(), syscall.SIGINT, syscall.SIGTERM)
	defer cancel()

	if config.timeout != 0 {
		ctx, cancel = context.WithTimeout(ctx, config.timeout)
		defer cancel()
	}

	test := testBase64Prefix(config.prefix)

	n, attempts, ok := findPointParallel(ctx, runtime.NumCPU(), p0, test, 0)

	private := "-"
	public := config.prefix + "..."
	if ok {
		scalarN := scalarFromUint64(n)

		po := new(edwards25519.Point).ScalarMult(scalarN, pointOffset)
		p := new(edwards25519.Point).Add(p0, po)
		public = base64.StdEncoding.EncodeToString(p.BytesMontgomery())

		if s0 != nil {
			so := new(edwards25519.Scalar).Multiply(scalarN, scalarOffset)
			s := new(edwards25519.Scalar).Add(s0, so)
			private = base64.StdEncoding.EncodeToString(scalarToKeyBytes(s))
		}
	}

	if config.output == "offset" {
		fmt.Printf("%d\n", n)
	} else {
		duration := time.Since(start)

		fmt.Printf("%-44s %-44s %-10s %-10s %s\n", "private", "public", "attempts", "duration", "attempts/s")
		fmt.Printf("%-44s %-44s %-10d %-10s %.0f\n", private, public, attempts, duration.Round(time.Second), float64(attempts)/duration.Seconds())
	}

	if !ok {
		os.Exit(1)
	}
}

func cmdAdd(args []string) {
	config := struct {
		offset uint64
		prefix string
	}{}

	fs := flag.NewFlagSet("add", flag.ExitOnError)
	fs.Uint64Var(&config.offset, "offset", 0, "add specified offset to the private key")
	fs.StringVar(&config.prefix, "prefix", "", "prefix of base64-encoded public key")
	fs.Parse(args)

	in, err := io.ReadAll(io.LimitReader(os.Stdin, 44))
	if err != nil {
		panic(err)
	}
	private := make([]byte, base64.StdEncoding.DecodedLen(len(in)))
	if n, err := base64.StdEncoding.Decode(private, in); err != nil {
		panic(err)
	} else if n != 32 {
		panic(fmt.Sprintf("Wrong private key length: %d", n))
	} else {
		private = private[:n]
	}

	s, _ := newPairFrom(private)

	// Worker starts search from a public key which corresponds to two points P and -P
	// but [parsePublicKey] returns (arbitrary?) one.
	// To find the right private key calculate both
	//     s + offset*scalarOffset
	//     s - offset*scalarOffset
	// and pick the one whose public key has the right prefix.
	so := edwards25519.NewScalar().Multiply(scalarFromUint64(config.offset), scalarOffset)

	sPlus := edwards25519.NewScalar().Add(s, so)
	sMinus := edwards25519.NewScalar().Subtract(s, so)

	skPlus := base64.StdEncoding.EncodeToString(scalarToKeyBytes(sPlus))
	skMinus := base64.StdEncoding.EncodeToString(scalarToKeyBytes(sMinus))

	pPlus := new(edwards25519.Point).ScalarBaseMult(sPlus)
	pkPlus := base64.StdEncoding.EncodeToString(pPlus.BytesMontgomery())

	pMinus := new(edwards25519.Point).ScalarBaseMult(sMinus)
	pkMinus := base64.StdEncoding.EncodeToString(pMinus.BytesMontgomery())

	if config.prefix != "" {
		if strings.HasPrefix(pkPlus, config.prefix) {
			fmt.Println(skPlus)
		} else if strings.HasPrefix(pkMinus, config.prefix) {
			fmt.Println(skMinus)
		} else {
			panic("prefix mismatch")
		}
	} else {
		fmt.Println(skPlus)
		fmt.Println(pkPlus)
		fmt.Println(skMinus)
		fmt.Println(pkMinus)
		os.Exit(1)
	}
}

func newPair() (*edwards25519.Scalar, *edwards25519.Point) {
	var key [32]byte
	if _, err := io.ReadFull(rand.Reader, key[:]); err != nil {
		panic(err)
	}
	return newPairFrom(key[:])
}

func newPairFrom(key []byte) (*edwards25519.Scalar, *edwards25519.Point) {
	s, err := edwards25519.NewScalar().SetBytesWithClamping(key)
	if err != nil {
		panic(err)
	}

	p := new(edwards25519.Point).ScalarBaseMult(s)

	return s, p
}

// parsePublicKey decodes base64-encoded public key
// and returns corresponding [edwards25519.Point] or error.
//
// https://datatracker.ietf.org/doc/html/rfc7748#section-4.1
func parsePublicKey(pk string) (*edwards25519.Point, error) {
	pkb, err := base64.StdEncoding.DecodeString(pk)
	if err != nil {
		return nil, err
	} else if len(pkb) != 32 {
		return nil, fmt.Errorf("wrong public key length")
	}

	u, _ := new(field.Element).SetBytes(pkb)

	// y = (u - 1) / (u + 1)
	var y, n, d, r field.Element

	n.Subtract(u, new(field.Element).One())
	d.Add(u, new(field.Element).One())
	r.Invert(&d)
	y.Multiply(&n, &r)

	return new(edwards25519.Point).SetBytes(y.Bytes())
}

func findPointParallel(ctx context.Context, workers int, p0 *edwards25519.Point, test func([]byte) bool, limit uint64) (uint64, uint64, bool) {
	result := make(chan uint64, workers)
	var attempts atomic.Uint64

	gctx, cancel := context.WithCancel(ctx)
	defer cancel()

	var wg sync.WaitGroup
	wg.Add(workers)
	for range workers {
		go func() {
			defer wg.Done()

			skip := randUint64()
			n, ok := findBatchPoint(gctx, p0, skip, 1024, test, limit/uint64(workers))

			attempts.Add(n - skip)
			if ok {
				result <- n
				cancel()
			}
		}()
	}
	wg.Wait()

	select {
	case n := <-result:
		return n, attempts.Load(), true
	case <-ctx.Done():
		return 0, attempts.Load(), false
	}
}

func findBatchPoint(ctx context.Context, p0 *edwards25519.Point, skip uint64, batchSize int, test func([]byte) bool, limit uint64) (uint64, bool) {
	skipOffset := new(edwards25519.Point).ScalarMult(scalarFromUint64(skip), pointOffset)
	p := new(edwards25519.Point).Add(p0, skipOffset)

	n := skip

	ua := make([]field.Element, batchSize+1)
	ub := make([]field.Element, batchSize+1)
	u := make([]field.Element, batchSize+1)
	scratch := make([][]field.Element, 2)

	for i := range scratch {
		scratch[i] = make([]field.Element, batchSize+1)
	}

	// offsets[i] = i * pointOffset
	offsets := make([]affine, batchSize)
	offsets[0].zero()
	poi := new(edwards25519.Point).Set(pointOffset)
	for i := 1; i < len(offsets); i++ {
		offsets[i].fromP3(poi)
		poi.Add(poi, pointOffset)
	}
	batchOffset := poi // batchSize * pointOffset

	// Starting point for the current batch
	pa := new(affine).fromP3(p)

	var bm [32]byte
	// Loop body complexity: (6M + 4A)*batchSize + 278M + 9A
	for {
		select {
		case <-ctx.Done():
			return n, false
		default:
		}

		// Complexity: (2M + 4A)*batchSize
		for i := range batchSize {
			sumU(pa, &offsets[i], &ua[i], &ub[i])
		}

		// Complexity: 9M + 9A
		p.Add(p, batchOffset)

		// Piggy-back on vector division to calculate p.zInv
		_, _, Z, _ := p.ExtendedCoordinates()
		ua[batchSize].One()
		ub[batchSize].Set(Z)

		// Complexity: 262M + 4M*(batchSize+1)
		vectorDivision(ua, ub, u, scratch)

		for i := range batchSize {
			copy(bm[:], u[i].Bytes()) // eliminate field.Element.Bytes() allocations
			if test(bm[:]) {
				return n + uint64(i), true
			}
		}

		n += uint64(batchSize)
		// Complexity: 3M
		pa.fromP3zInv(p, &u[batchSize])

		if limit > 0 {
			if limit <= uint64(batchSize) {
				// TODO: should signal finish somehow else
				return n, true
			}
			limit -= uint64(batchSize)
		}
	}
}

// sumU returns u-Montgomery coordinate of p + q as a fraction u = ua / ub
// Complexity: 2M + 4A
func sumU(p, q *affine, ua, ub *field.Element) {
	// Affine addition formulae (independent of d) for twisted Edwards curves
	// https://eprint.iacr.org/2008/522.pdf
	//
	// y3 = (x1*y1 − x2*y2) / (x1*y2 − y1*x2) = nom/den
	var nom, den field.Element

	nom.Subtract(&p.XY, &q.XY)
	den.Subtract(
		new(field.Element).Multiply(&p.X, &q.Y),
		new(field.Element).Multiply(&p.Y, &q.X),
	)

	// https://www.rfc-editor.org/rfc/rfc7748.html#section-4.1
	// (u, v) = ((1+y)/(1-y), sqrt(-486664)*u/x)
	// u = (1 + y) / (1 - y) = (1 + nom/den) / (1 - nom/den) = (den + nom) / (den - nom) = ua/ub
	ua.Add(&den, &nom)
	ub.Subtract(&den, &nom)
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

func testBase64Prefix(prefix string) func([]byte) bool {
	decoded, decodedBits := decodeBase64PrefixBits(prefix)

	if decodedBits%8 == 0 {
		return func(b []byte) bool {
			return bytes.HasPrefix(b, decoded)
		}
	}

	decodedBytes := decodedBits / 8
	shift := 8 - (decodedBits % 8)
	tailByte := decoded[decodedBytes] >> shift
	decoded = decoded[:decodedBytes]

	return func(b []byte) bool {
		return len(b) > decodedBytes && // must be long enough to check tail byte
			bytes.Equal(b[:decodedBytes], decoded) &&
			b[decodedBytes]>>shift == tailByte
	}
}

// decodeBase64PrefixBits returns decoded prefix and number of decoded bits.
func decodeBase64PrefixBits(prefix string) ([]byte, int) {
	decodedBits := 6 * len(prefix)
	tailBits := decodedBits % 8

	// Parse prefix as base64 number
	decodedInt := big.NewInt(0)
	_64 := big.NewInt(64)
	const stdEncoding = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/"
	for i := 0; i < len(prefix); i++ {
		if digit := strings.IndexByte(stdEncoding, prefix[i]); digit >= 0 {
			decodedInt.Mul(decodedInt, _64)
			decodedInt.Add(decodedInt, big.NewInt(int64(digit)))
		} else if prefix[i] == '=' {
			// skip
		} else {
			panic("invalid base64 byte: " + prefix[i:i+1])
		}
	}
	if tailBits != 0 {
		decodedInt.Mul(decodedInt, big.NewInt(1<<(8-tailBits)))
	}
	decoded := decodedInt.Bytes()

	// left pad decoded prefix with zeros
	buf := make([]byte, (decodedBits+7)/8)
	copy(buf[len(buf)-len(decoded):], decoded)

	return buf, decodedBits
}

func randUint64() uint64 {
	var num uint64
	err := binary.Read(rand.Reader, binary.NativeEndian, &num)
	if err != nil {
		panic(err)
	}
	return num
}

// vectorDivision calculates u = x / y using scratch.
//
// vectorDivision uses:
//
//	4*(n-1)+1 multiplications
//	1 invert = ~265 multiplications
//
// Complexity: 262M + 4M*n
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

type affine struct {
	X, Y, XY field.Element
}

func (v *affine) zero() *affine {
	v.X.Zero()
	v.Y.One()
	v.XY.Zero()
	return v
}

// Complexity: 1I + 3M = 268M
func (v *affine) fromP3(p *edwards25519.Point) *affine {
	X, Y, Z, _ := p.ExtendedCoordinates()
	var zInv field.Element
	zInv.Invert(Z)
	v.X.Multiply(X, &zInv)
	v.Y.Multiply(Y, &zInv)
	v.XY.Multiply(&v.X, &v.Y)
	return v
}

// Complexity: 3M
func (v *affine) fromP3zInv(p *edwards25519.Point, zInv *field.Element) *affine {
	X, Y, _, _ := p.ExtendedCoordinates()
	v.X.Multiply(X, zInv)
	v.Y.Multiply(Y, zInv)
	v.XY.Multiply(&v.X, &v.Y)
	return v
}

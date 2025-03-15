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

// Ed25519 group's cofactor
const ed25519GroupCofactor = 8

var (
	scalarOffset = scalarFromBytes(ed25519GroupCofactor)
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

	var s0b []byte
	var p0 *edwards25519.Point
	var err error

	if config.public != "" {
		p0, err = parsePublicKey(config.public)
		if err != nil {
			panic(err)
		}
	} else {
		s0b = make([]byte, 32)
		if _, err := io.ReadFull(rand.Reader, s0b); err != nil {
			panic(err)
		}
		clampKeyBytes(s0b)
		_, p0 = newPairFrom(s0b)
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

		if len(s0b) > 0 {
			s0 := fieldElementFromBytes(s0b)
			so := fieldElementFromUint64(n)
			so.Mult32(so, ed25519GroupCofactor)
			s := s0.Add(s0, so)

			private = base64.StdEncoding.EncodeToString(s.Bytes())
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

	in := make([]byte, 44)
	if _, err := io.ReadFull(os.Stdin, in); err != nil {
		panic(err)
	}

	s0b := make([]byte, 32)
	if n, err := base64.StdEncoding.Decode(s0b, in); err != nil {
		panic(err)
	} else if n != 32 {
		panic(fmt.Sprintf("Wrong private key length: %d", n))
	}
	clampKeyBytes(s0b)

	s0 := fieldElementFromBytes(s0b)
	so := fieldElementFromUint64(config.offset)
	so.Mult32(so, ed25519GroupCofactor)

	// Worker starts search from a public key which corresponds to two points P and -P
	// but [parsePublicKey] returns (arbitrary?) one.
	// To find the right private key calculate both:
	//     s0 + offset*scalarOffset
	//     s0 - offset*scalarOffset
	// and pick the one whose public key has the right prefix.
	sPlus := new(field.Element).Add(s0, so)
	sMinus := new(field.Element).Subtract(s0, so)

	_, pPlus := newPairFrom(sPlus.Bytes())
	pkPlus := base64.StdEncoding.EncodeToString(pPlus.BytesMontgomery())
	_, pMinus := newPairFrom(sMinus.Bytes())
	pkMinus := base64.StdEncoding.EncodeToString(pMinus.BytesMontgomery())

	skPlus := base64.StdEncoding.EncodeToString(sPlus.Bytes())
	skMinus := base64.StdEncoding.EncodeToString(sMinus.Bytes())

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
	return s, new(edwards25519.Point).ScalarBaseMult(s)
}

// clampKeyBytes applies the buffer pruning described in RFC 8032,
// Section 5.1.5 (also known as clamping) and returns true if key has changed.
func clampKeyBytes(key []byte) bool {
	k0, k31 := key[0], key[31]
	key[0] &= 248
	key[31] &= 63
	key[31] |= 64
	return k0 != key[0] || k31 != key[31]
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

	u, err := new(field.Element).SetBytes(pkb)
	if err != nil {
		return nil, err
	}

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
	if batchSize < 0 || batchSize%2 != 0 {
		panic("batchSize must be non-negative and even")
	}

	skipOffset := new(edwards25519.Point).ScalarMult(scalarFromUint64(skip), pointOffset)
	p := new(edwards25519.Point).Add(p0, skipOffset)

	n := skip

	ua := make([]field.Element, batchSize+2)
	ub := make([]field.Element, batchSize+2)
	u := make([]field.Element, batchSize+2)
	scratch := make([][]field.Element, 2)

	for i := range scratch {
		scratch[i] = make([]field.Element, batchSize+2)
	}

	// offsets[i] = (i+1) * pointOffset
	offsets := make([]affine, batchSize/2)
	poi := new(edwards25519.Point).Set(pointOffset)
	for i := range batchSize/2 - 1 {
		offsets[i].fromP3(poi)
		poi.Add(poi, pointOffset)
	}
	// batchOffset = (batchSize+1) * pointOffset
	batchOffset := new(edwards25519.Point).Set(pointOffset)
	if batchSize > 0 {
		offsets[batchSize/2-1].fromP3(poi)
		batchOffset.Add(batchOffset, poi)
		batchOffset.Add(batchOffset, poi)
		p.Add(p, poi)
	}

	// Center point for the current batch
	pa := new(affine).fromP3(p)

	var bm [32]byte
	// One iteration tests batchSize+1 points of the
	// batch = {p − batchSize/2*pointOffset, ... , p − pointOffset, p, p + pointOffset, ... , p + batchSize/2*pointOffset}
	//
	// Complexity: (5M + 4A)*batchSize + 282M + 11A
	for {
		select {
		case <-ctx.Done():
			return n, false
		default:
		}

		// Calculate u-coordinates of (pa + offsets[i]), (pa − offsets[i]) and pa
		// points on a Montgomery curve as a ratio u = ua / ub
		//
		// Affine addition formulae (independent of d) for twisted Edwards curves,
		// see https://eprint.iacr.org/2008/522.pdf
		//
		// y3 = (x1*y1 − x2*y2) / (x1*y2 − y1*x2) = nom / den
		//
		// Symmetric negative point p2' = −p2 has y2' = y2 and x2' = −x2, therefore
		//
		// y3' = (x1*y1 + x2*y2) / (x1*y2 + y1*x2)
		//
		// The u-coordinate on a Montgomery curve,
		// see https://www.rfc-editor.org/rfc/rfc7748.html#section-4.1
		//
		// u = (1 + y) / (1 − y) = (1 + nom/den) / (1 − nom/den) = (den + nom) / (den − nom) = ua / ub
		//
		// Complexity: (2M + 8A)*batchSize/2 + 2A = (1M + 4A)*batchSize + 2A
		for i := range batchSize / 2 {
			p1, p2 := pa, &offsets[i]

			var nom, den, x1y2, y1x2 field.Element

			x1y2.Multiply(&p1.X, &p2.Y)
			y1x2.Multiply(&p1.Y, &p2.X)

			// p3 = p1 + p2
			// y3 = (x1*y1 − x2*y2) / (x1*y2 − y1*x2) = nom / den
			// u = (den + nom) / (den − nom)
			nom.Subtract(&p1.XY, &p2.XY)
			den.Subtract(&x1y2, &y1x2)
			ua[batchSize/2+1+i].Add(&den, &nom)
			ub[batchSize/2+1+i].Subtract(&den, &nom)

			// p3' = p1 − p2
			// y3' = (x1*y1 + x2*y2) / (x1*y2 + y1*x2) = nom / den
			// u' = (den + nom) / (den − nom)
			nom.Add(&p1.XY, &p2.XY)
			den.Add(&x1y2, &y1x2)
			ua[batchSize/2-1-i].Add(&den, &nom)
			ub[batchSize/2-1-i].Subtract(&den, &nom)
		}
		// pa is the center point of the batch
		ua[batchSize/2].Add(new(field.Element).One(), &pa.Y)
		ub[batchSize/2].Subtract(new(field.Element).One(), &pa.Y)

		// Complexity: 9M + 9A
		p.Add(p, batchOffset)

		// Piggyback on vector division to calculate 1/p.Z
		_, _, pZ, _ := p.ExtendedCoordinates()
		ua[batchSize+1].One()
		ub[batchSize+1].Set(pZ)
		pZinv := &u[batchSize+1]

		// Complexity: 262M + 4M*(batchSize+2) = 4M*batchSize + 270M
		vectorDivision(ua, ub, u, scratch)

		for i := range batchSize + 1 {
			copy(bm[:], u[i].Bytes()) // eliminate field.Element.Bytes() allocations
			if test(bm[:]) {
				return n + uint64(i), true
			}
		}

		n += uint64(batchSize + 1)
		// Complexity: 3M
		pa.fromP3zInv(p, pZinv)

		if limit > 0 {
			if limit <= uint64(batchSize) {
				// TODO: should signal finish somehow else
				return n, true
			}
			limit -= uint64(batchSize)
		}
	}
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

func fieldElementFromBytes(x []byte) *field.Element {
	var buf [32]byte
	copy(buf[:], x)
	fe, err := new(field.Element).SetBytes(buf[:])
	if err != nil {
		panic(err)
	}
	return fe
}

func fieldElementFromUint64(n uint64) *field.Element {
	var nb [8]byte
	binary.LittleEndian.PutUint64(nb[:], n)
	return fieldElementFromBytes(nb[:])
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

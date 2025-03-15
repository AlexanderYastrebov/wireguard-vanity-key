package main

import (
	"context"
	"crypto/ecdh"
	"crypto/rand"
	"encoding/base64"
	"fmt"
	"io"
	"math"
	"reflect"
	"runtime"
	"strings"
	"testing"

	"filippo.io/edwards25519"
)

func assertEqual(tb testing.TB, expected, got any) {
	tb.Helper()
	if !reflect.DeepEqual(expected, got) {
		tb.Errorf("%v != %v", expected, got)
	}
}

func requireEqual(tb testing.TB, expected, got any) {
	tb.Helper()
	if !reflect.DeepEqual(expected, got) {
		tb.Fatalf("%v != %v", expected, got)
	}
}

func BenchmarkNewPrivateKey(b *testing.B) {
	var key [32]byte
	_, err := io.ReadFull(rand.Reader, key[:])
	requireEqual(b, nil, err)

	b.ResetTimer()
	for range b.N {
		priv, _ := ecdh.X25519().NewPrivateKey(key[:])
		_ = priv.PublicKey().Bytes()
	}
}

var testPrefix = testBase64Prefix("GoodLuckWithThisPrefix")

func BenchmarkFindBatchPoint(b *testing.B) {
	for _, batchSize := range []int{
		2, 32, 64, 128, 256, 512, 1024,
		2048, 4096, 8192,
	} {
		b.Run(fmt.Sprintf("%d", batchSize), func(b *testing.B) {
			_, p0 := newPair()

			ctx, cancel := context.WithCancel(context.Background())
			defer cancel()

			skip := randUint64()

			b.ResetTimer()
			findBatchPoint(ctx, p0, skip, batchSize, testPrefix, uint64(b.N))
		})
	}
}

func BenchmarkFindPointParallel(b *testing.B) {
	_, p0 := newPair()

	b.ResetTimer()
	findPointParallel(context.Background(), min(runtime.NumCPU(), b.N), p0, testPrefix, uint64(b.N))
}

func TestTestBase64Prefix(t *testing.T) {
	plain := []byte("Hello World!")
	encoded := base64.StdEncoding.EncodeToString([]byte(plain))

	for i := 1; i <= len(encoded); i++ {
		prefix := encoded[:i]

		t.Run(prefix, func(t *testing.T) {
			test := testBase64Prefix(prefix)
			if !test([]byte(plain)) {
				t.Errorf("Prefix mismatch: %q %q", plain, prefix)
			}
		})
	}

	// Single symbol is a 6-bit prefix
	assertEqual(t, true, testBase64Prefix("A")([]byte{0}))
	assertEqual(t, true, testBase64Prefix("A")([]byte{0b000000_01}))
	assertEqual(t, true, testBase64Prefix("A")([]byte{0b000000_10}))
	assertEqual(t, true, testBase64Prefix("A")([]byte{0b000000_11}))
	assertEqual(t, false, testBase64Prefix("A")([]byte{0b000001_00}))

	assertEqual(t, true, testBase64Prefix("B")([]byte{0b000001_00}))
	assertEqual(t, true, testBase64Prefix("B")([]byte{0b000001_01}))
	assertEqual(t, true, testBase64Prefix("B")([]byte{0b000001_10}))
	assertEqual(t, true, testBase64Prefix("B")([]byte{0b000001_11}))
	assertEqual(t, false, testBase64Prefix("B")([]byte{0}))
	assertEqual(t, false, testBase64Prefix("B")([]byte{1}))

	// Two symbols is a 12-bit prefix
	assertEqual(t, true, testBase64Prefix("AA")([]byte{0, 0}))
	assertEqual(t, true, testBase64Prefix("AA")([]byte{0, 0b0000_0001}))
	assertEqual(t, true, testBase64Prefix("AB")([]byte{0b000000_00, 0b0001_0000}))
	assertEqual(t, true, testBase64Prefix("AB")([]byte{0b000000_00, 0b0001_0001}))
	assertEqual(t, true, testBase64Prefix("AB")([]byte{0b000000_00, 0b0001_0010}))
	assertEqual(t, true, testBase64Prefix("BB")([]byte{0b000001_00, 0b0001_0000}))
	assertEqual(t, false, testBase64Prefix("BB")([]byte{0b000001_01, 0b0001_0000}))

	assertEqual(t, true, testBase64Prefix("AAA")([]byte{0, 0, 0}))
	assertEqual(t, true, testBase64Prefix("AAA")([]byte{0, 0, 0b00_000001}))
}

func TestParsePublicKey(t *testing.T) {
	for _, pk := range []string{
		"QiyOemIn17yhNQs+K7cnn3iXuHu2hUt4PGDoAxGuMHk=",
		"vQnB//PF0URzwwsH0b1ff7a0P3jLKbrOCdLiTkWkvQA=",
		"3nN+Tj4J/e99YWD6TFMvhfMNJCrORoSf8ommtXeXvBs=",
		"Fo8iOSvqtfDjtBALpwGALNiwaZNgMrQYXIEDB2oU6lQ=",
		"YR3nSufwy4r5FuCE7GujLSLssyVJ6iKy2utbUCQelh4=",
	} {
		t.Run(pk, func(t *testing.T) {
			p, err := parsePublicKey(pk)
			requireEqual(t, nil, err)

			assertEqual(t, pk, base64.StdEncoding.EncodeToString(p.BytesMontgomery()))
		})
	}
}

func TestFindBatchPoint(t *testing.T) {
	t.Run("qkHBetbXfAxsmr0jH6Zs6Dx1ZEReO9WBZCoNREce0gE=", func(t *testing.T) {
		p0, err := parsePublicKey("qkHBetbXfAxsmr0jH6Zs6Dx1ZEReO9WBZCoNREce0gE=")
		requireEqual(t, nil, err)

		const expectedOffset uint64 = 92950

		offset, ok := findBatchPoint(context.Background(), p0, 0, 1024, testBase64Prefix("AY/"), 0)
		assertEqual(t, true, ok)
		assertEqual(t, expectedOffset, offset)

		p := new(edwards25519.Point).Add(p0, new(edwards25519.Point).ScalarMult(scalarFromUint64(offset), pointOffset))
		assertEqual(t, "AY/yq7zukqRmMUzqqPFmtqXJdAcbmh8mn4rMgtjVnGI=", base64.StdEncoding.EncodeToString(p.BytesMontgomery()))
	})

	t.Run("params", func(t *testing.T) {
		// TODO: fix overflow skip+offset
		for _, skip := range []uint64{0, 1, 2, 3, 10, math.MaxUint64 / 2} {
			for _, batchSize := range []int{0, 2, 512, 1024} {
				t.Run(fmt.Sprintf("skip=%d,batchSize=%d", skip, batchSize), func(t *testing.T) {
					s0, p0 := newPair()
					t.Logf("s0: %s", base64.StdEncoding.EncodeToString(s0.Bytes()))
					t.Logf("p0: %s", base64.StdEncoding.EncodeToString(p0.BytesMontgomery()))

					offset, ok := findBatchPoint(context.Background(), p0, skip, batchSize, testBase64Prefix("AY/"), 0)
					assertEqual(t, true, ok)

					so := new(edwards25519.Scalar).Multiply(scalarFromUint64(offset), scalarOffset)
					s := new(edwards25519.Scalar).Add(s0, so)
					p := new(edwards25519.Point).ScalarBaseMult(s)
					t.Logf("s: %s", base64.StdEncoding.EncodeToString(s0.Bytes()))
					t.Logf("p: %s", base64.StdEncoding.EncodeToString(p0.BytesMontgomery()))

					assertEqual(t, true, strings.HasPrefix(base64.StdEncoding.EncodeToString(p.BytesMontgomery()), "AY/"))
				})
			}
		}
	})
}

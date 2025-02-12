package main

import (
	"bytes"
	"context"
	"crypto/ecdh"
	"crypto/rand"
	"encoding/base64"
	"fmt"
	"io"
	"runtime"
	"sync/atomic"
	"testing"

	"filippo.io/edwards25519"
	"filippo.io/edwards25519/field"
)

func BenchmarkNewPrivateKey(b *testing.B) {
	var key [32]byte
	_, err := io.ReadFull(rand.Reader, key[:])
	if err != nil {
		b.Fatal(err)
	}
	b.ResetTimer()

	for range b.N {
		priv, _ := ecdh.X25519().NewPrivateKey(key[:])
		_ = priv.PublicKey().Bytes()
	}
}

var testPrefix = testBase64Prefix("GoodLuckWithThisPrefix")

func BenchmarkFindBatchPoint(b *testing.B) {
	for _, batchSize := range []int{
		1, 32, 64, 128, 256, 512, 1024,
		2048, 4096, 8192,
	} {
		b.Run(fmt.Sprintf("%d", batchSize), func(b *testing.B) {
			_, p0 := newPair()

			i := b.N

			ctx, cancel := context.WithCancel(context.Background())
			defer cancel()

			b.ResetTimer()
			findBatchPoint(ctx, p0, randUint64(), batchSize, func(p []byte) bool {
				match := testPrefix(p)
				i--
				return i == 0 || match
			})
		})
	}
}

func BenchmarkFindPointParallel(b *testing.B) {
	_, p0 := newPair()

	var i atomic.Int64
	i.Store(int64(b.N))

	findPointParallel(context.Background(), runtime.NumCPU(), p0, func(p []byte) bool {
		match := testPrefix(p)
		return i.Add(-1) <= 0 || match
	})
}

func TestBatchBytesMontgomery(t *testing.T) {
	pts := make([]edwards25519.Point, 64)
	u := make([]field.Element, len(pts))
	scratch := make([][]field.Element, 4)

	for i := range scratch {
		scratch[i] = make([]field.Element, len(pts))
	}

	for i := range pts {
		_, p := newPair()
		pts[i].Set(p)
	}

	batchBytesMontgomery(pts, u, scratch)

	for i, p := range pts {
		if !bytes.Equal(p.BytesMontgomery(), u[i].Bytes()) {
			t.Errorf("Wrong montgomery bytes")
		}
	}

	t.Run("no allocs", func(t *testing.T) {
		n := testing.AllocsPerRun(100, func() {
			batchBytesMontgomery(pts, u, scratch)
		})
		if n != 0 {
			t.Errorf("Unexpected allocations: %.0f", n)
		}
	})
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

	assertEqual := func(a, b bool) {
		t.Helper()
		if a != b {
			t.Error("<-- see")
		}
	}

	// Single symbol is a 6-bit prefix
	assertEqual(true, testBase64Prefix("A")([]byte{0}))
	assertEqual(true, testBase64Prefix("A")([]byte{0b000000_01}))
	assertEqual(true, testBase64Prefix("A")([]byte{0b000000_10}))
	assertEqual(true, testBase64Prefix("A")([]byte{0b000000_11}))
	assertEqual(false, testBase64Prefix("A")([]byte{0b000001_00}))

	assertEqual(true, testBase64Prefix("B")([]byte{0b000001_00}))
	assertEqual(true, testBase64Prefix("B")([]byte{0b000001_01}))
	assertEqual(true, testBase64Prefix("B")([]byte{0b000001_10}))
	assertEqual(true, testBase64Prefix("B")([]byte{0b000001_11}))
	assertEqual(false, testBase64Prefix("B")([]byte{0}))
	assertEqual(false, testBase64Prefix("B")([]byte{1}))

	// Two symbols is a 12-bit prefix
	assertEqual(true, testBase64Prefix("AA")([]byte{0, 0}))
	assertEqual(true, testBase64Prefix("AA")([]byte{0, 0b0000_0001}))
	assertEqual(true, testBase64Prefix("AB")([]byte{0b000000_00, 0b0001_0000}))
	assertEqual(true, testBase64Prefix("AB")([]byte{0b000000_00, 0b0001_0001}))
	assertEqual(true, testBase64Prefix("AB")([]byte{0b000000_00, 0b0001_0010}))
	assertEqual(true, testBase64Prefix("BB")([]byte{0b000001_00, 0b0001_0000}))
	assertEqual(false, testBase64Prefix("BB")([]byte{0b000001_01, 0b0001_0000}))

	assertEqual(true, testBase64Prefix("AAA")([]byte{0, 0, 0}))
	assertEqual(true, testBase64Prefix("AAA")([]byte{0, 0, 0b00_000001}))
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
			if err != nil {
				t.Fatal(err)
			}

			if got := base64.StdEncoding.EncodeToString(p.BytesMontgomery()); got != pk {
				t.Errorf("pk: %s, got: %s", pk, got)
			}
		})
	}
}

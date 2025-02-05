package main

import (
	"bytes"
	"context"
	"crypto/ecdh"
	"crypto/rand"
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

func BenchmarkFindPoint(b *testing.B) {
	_, p0 := newPair()

	i := b.N

	findPoint(context.Background(), p0, randUint64(), func(p []byte) bool {
		match := hasBase64Prefix(p, []byte("GoodLuckWithThisPrefix"))
		i--
		return i == 0 || match
	})
}

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
				match := hasBase64Prefix(p, []byte("GoodLuckWithThisPrefix"))
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
		match := hasBase64Prefix(p, []byte("GoodLuckWithThisPrefix"))
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

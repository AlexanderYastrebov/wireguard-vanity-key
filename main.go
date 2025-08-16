// Package main searches for a [WireGuard] [curve25519] keypair
// with a base64-encoded public key that has a specified prefix.
//
// [WireGuard]: https://www.wireguard.com/
// [curve25519]: https://datatracker.ietf.org/doc/html/rfc7748#section-4.1
package main

import (
	"context"
	"crypto/ecdh"
	"crypto/rand"
	"encoding/base64"
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

	"github.com/AlexanderYastrebov/vanity25519"
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

	var startKey *ecdh.PrivateKey
	var startPublicKey []byte
	var err error

	if config.public != "" {
		startPublicKey, err = base64.StdEncoding.DecodeString(config.public)
		if err != nil {
			panic(err)
		}
	} else {
		startKey, err = ecdh.X25519().GenerateKey(rand.Reader)
		if err != nil {
			panic(err)
		}
		startPublicKey = startKey.PublicKey().Bytes()
	}

	ctx, cancel := signal.NotifyContext(context.Background(), syscall.SIGINT, syscall.SIGTERM)
	defer cancel()

	if config.timeout != 0 {
		ctx, cancel = context.WithTimeout(ctx, config.timeout)
		defer cancel()
	}

	test := vanity25519.HasPrefixBits(decodeBase64PrefixBits(config.prefix))

	vanityPublicKey, offset, attempts, ok := searchParallel(ctx, runtime.NumCPU(), startPublicKey, test)

	private := "-"
	public := config.prefix + "..."
	if ok {
		public = base64.StdEncoding.EncodeToString(vanityPublicKey)
		if startKey != nil {
			vanityPrivateKey, err := vanity25519.Add(startKey.Bytes(), offset)
			if err != nil {
				panic(err)
			}
			private = base64.StdEncoding.EncodeToString(vanityPrivateKey)
		}
	}

	if config.output == "offset" {
		fmt.Println(offset.String())
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
		offset *big.Int
	}{}
	var ok bool

	fs := flag.NewFlagSet("add", flag.ExitOnError)
	fs.Func("offset", "add specified offset to the private key", func(s string) error {
		if config.offset, ok = new(big.Int).SetString(s, 10); !ok {
			return fmt.Errorf("invalid offset")
		}
		return nil
	})
	fs.Parse(args)

	if config.offset == nil {
		panic("offset required")
	}

	in := make([]byte, 44)
	if _, err := io.ReadFull(os.Stdin, in); err != nil {
		panic(err)
	}
	startPrivateKey, err := base64.StdEncoding.DecodeString(string(in))
	if err != nil {
		panic(err)
	}

	vanityPrivateKey, err := vanity25519.Add(startPrivateKey, config.offset)
	if err != nil {
		panic(err)
	}
	fmt.Println(base64.StdEncoding.EncodeToString(vanityPrivateKey))
}

func searchParallel(ctx context.Context, workers int, startPublicKey []byte, test func([]byte) bool) ([]byte, *big.Int, uint64, bool) {
	type result struct {
		publicKey []byte
		offset    *big.Int
	}
	found := make(chan result, workers)
	var totalAttempts atomic.Uint64

	gtx, cancel := context.WithCancel(ctx)
	defer cancel()

	var wg sync.WaitGroup
	for range workers {
		wg.Go(func() {
			var attempts uint64
			testAttempts := func(b []byte) bool {
				attempts++
				return test(b)
			}

			vanity25519.Search(gtx, startPublicKey, randBigInt(), 4096, testAttempts, func(publicKey []byte, offset *big.Int) {
				found <- result{publicKey, offset}
				cancel()
			})

			totalAttempts.Add(attempts)
		})
	}
	wg.Wait()

	select {
	case r := <-found:
		return r.publicKey, r.offset, totalAttempts.Load(), true
	case <-ctx.Done():
		return nil, big.NewInt(0), totalAttempts.Load(), false
	}
}

// decodeBase64PrefixBits returns decoded prefix and number of decoded bits.
func decodeBase64PrefixBits(prefix string) ([]byte, int) {
	decodedBits := 6 * len(prefix)
	quantums := (len(prefix) + 3) / 4
	prefix += strings.Repeat("A", quantums*4-len(prefix))
	buf := make([]byte, quantums*3)
	_, err := base64.StdEncoding.Decode(buf, []byte(prefix))
	if err != nil {
		panic(err)
	}
	return buf, decodedBits
}

func randBigInt() *big.Int {
	var buf [8]byte
	rand.Read(buf[:])
	return new(big.Int).SetBytes(buf[:])
}

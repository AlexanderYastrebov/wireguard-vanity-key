# wireguard-vanity-key

This tool searches for a [WireGuard](https://www.wireguard.com/) Curve25519 keypair
with a base64-encoded public key that has a specified prefix.

Compared to [similar tools](#similar-tools), it uses [the fastest search algorithm](#the-fastest-search-algorithm) 🚀

## Example

Install the tool locally and run:
```console
$ go install github.com/AlexanderYastrebov/wireguard-vanity-key@latest
$ wireguard-vanity-key -prefix=2025
private                                      public                                       attempts   duration   attempts/s
WHakaGFouuy2AxMmOdSTf2L2KWsI6a3s+gvAOKuKtH0= 2025sb38RUVI+GJg5Uk2RRPuJfhZyg4uSxfV2WDn1g8= 47423039   2s         19926032

$ # verify
$ echo WHakaGFouuy2AxMmOdSTf2L2KWsI6a3s+gvAOKuKtH0= | wg pubkey
2025sb38RUVI+GJg5Uk2RRPuJfhZyg4uSxfV2WDn1g8=
```

or run the tool from the source repository:
```console
$ go run . -prefix=2025
```

or use the Docker image:
```console
$ docker pull ghcr.io/alexanderyastrebov/wireguard-vanity-key:latest
$ docker run  ghcr.io/alexanderyastrebov/wireguard-vanity-key:latest -prefix=2025
```

## Performance

The tool checks ~29'000'000 keys per second on a test machine:

```console
$ go run . -prefix=GoodLuckWithThisPrefix -timeout=20s
private                                      public                                       attempts   duration   attempts/s
-                                            GoodLuckWithThisPrefix...                    583920640  20s        29194831
```

In practice, it finds a 4-character prefix in a second and a 5-character prefix in a minute:
```console
$ while go run . -prefix=AYAYA ; do : ; done
private                                      public                                       attempts   duration   attempts/s
OM9WvIxO90NRnHpMLBYbKCwRxj1KcwWVfo5EO1vftls= AYAYAgcnXbdsMwLB+nR0kkWDpIkMr+3thhfGnBEvTmM= 515859548  23s        22328225
private                                      public                                       attempts   duration   attempts/s
eEbiqUhcUrH6Uj1p7cycgTOspY6fMxxImSNNr1YvaEg= AYAYA4yow92Ks1wnbQeceKEWIYHhaRyezomUz9SQJic= 350598404  19s        18060407
private                                      public                                       attempts   duration   attempts/s
OOkjEu4elrWJ4MD+OxB2kvUcKdyo482E3G3Y/tLBsmI= AYAYAW4yGEUVT/IkX3T6ZZTnz3yPS1lPxiRe0yhOCCs= 260273230  14s        17972036
private                                      public                                       attempts   duration   attempts/s
+BWkcGvbkXFxNgxIrAYyJoMF1R6R3eguv5NyMsdlaEA= AYAYAQEsY0gagwZ5lGLRQYfxQ+5rl83LOPmaASvASFQ= 1094012149 56s        19446702
private                                      public                                       attempts   duration   attempts/s
aG7Rakjbn1kpc2HN7fUz1u/ZrTcYziXg7OJq2EcMWFU= AYAYAe9QZdXn36CrkOK8aoD8h92mbEHCQt1QdTBARjY= 1088287697 56s        19483959
```

Each additional character increases search time by a factor of 64.

## Blind search

The tool supports blind search, i.e., when the worker does not know the private key. See [demo-blind.sh](demo-blind.sh).

## Similar tools

* [wireguard-vanity-address](https://github.com/warner/wireguard-vanity-address)
* [wireguard-vanity-keygen](https://github.com/axllent/wireguard-vanity-keygen)
* [Wireguard-Vanity-Key-Searcher](https://github.com/volleybus/Wireguard-Vanity-Key-Searcher)
* [wgmine](https://github.com/thatsed/wgmine)
* [Vanity](https://github.com/samuel-lucas6/Vanity)

## The fastest search algorithm

A WireGuard key pair consists of a 256-bit random private key and a public key derived by scalar multiplication on [Curve25519](https://en.wikipedia.org/wiki/Curve25519) involving arithmetic operations (additions, multiplications) in a finite field.

The performance of any brute-force key search algorithm ultimately depends on the number of finite field **multiplications**
per candidate key - the most expensive field operation.

All available WireGuard vanity key search tools use the straightforward approach:
multiply the base point by a random candidate private key and check the resulting public key:
```
public_key = private_key × base_point
```
For the WireGuard key format, this basic algorithm requires **2561** field multiplications (using [square-and-multiply](https://github.com/golang/go/commit/e005cdc62081130117a3fa30d01cd28ee076ed93)) or **743** field multiplications (using [Twisted Edwards curve](https://github.com/FiloSottile/edwards25519/commit/2941d4c8cdacb392a1b39f85adafaeae65bb50f6)) per candidate key.

This tool uses only **5 (five)** field multiplications per candidate key and other optimizations, which makes it the fastest 🚀

### ➕ Point increment approach

Inspired by [wireguard-vanity-address "faster algorithm"](https://github.com/warner/wireguard-vanity-address/pull/15),
instead of doing full scalar multiplication for each candidate, this tool applies a point increment technique that reduces the number of multiplications:
```
public_key0 = private_key0 × base_point
public_key1 = (private_key0 + const_offset) × base_point 
            = private_key0 × base_point + const_offset × base_point 
            = public_key0 + const_offset × base_point
            = public_key0 + const_point_offset
```
i.e., the candidate public key is obtained by point addition instead of base point multiplication,
which requires fewer field multiplications (275 vs. 743 or ~60% faster).
See [initial commit](https://github.com/AlexanderYastrebov/wireguard-vanity-key/commit/8c25defadec12585a80245faa40ddda2c192d423).

### 🧮 Batch field inversion

Getting the public key in WireGuard format (Montgomery form) requires a field inversion, which uses **265** field multiplications.
This tool uses the Montgomery trick to implement batch inversion, so for a batch of N candidates,
only **1 (one)** inversion is needed, resulting in a huge speedup.
See https://github.com/AlexanderYastrebov/wireguard-vanity-key/pull/3.

### 🪞 Affine coordinates and offset symmetry

The algorithm uses affine coordinates and exploits symmetries in precomputed point offsets,
saving even more field multiplications. See https://github.com/AlexanderYastrebov/wireguard-vanity-key/pull/10, https://github.com/AlexanderYastrebov/wireguard-vanity-key/pull/12, and
https://github.com/AlexanderYastrebov/wireguard-vanity-key/pull/14.

### ⚡ Fast base64 prefix check

Other tools encode the full public key to base64 and compare the prefix. This tool decodes the base64 prefix and compares it to public key bytes directly. See https://github.com/AlexanderYastrebov/wireguard-vanity-key/pull/5.

### 🏆 High-performance C implementation

For raw speed, the C worker (`wvk`) uses [awslabs/s2n-bignum](https://github.com/awslabs/s2n-bignum) - 
a highly optimized field arithmetic library written in assembly.
The worker supports prefix lengths up to 10 base64 characters, so the prefix check becomes a single masked integer comparison.
These two optimizations make `wvk` ~2 times faster than the Go implementation.
See https://github.com/AlexanderYastrebov/wireguard-vanity-key/pull/15.

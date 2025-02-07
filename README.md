# wireguard-vanity-key

Inspired by [wireguard-vanity-address "faster algorithm"](https://github.com/warner/wireguard-vanity-address/pull/15),
this tool searches for a [WireGuard](https://www.wireguard.com/) Curve25519 keypair
with a base64-encoded public key that has a specified prefix.

## Example

Install the tool locally and run:
```console
$ go install github.com/AlexanderYastrebov/wireguard-vanity-key@latest
$ wireguard-vanity-key --prefix=2025
private                                      public                                       attempts   duration   attempts/s
4JFWFevraBg5yLD2rCFzGMpbKZOC3BV5rNom+Um7EGg= 2025Q6KcDb+v/nj2/ErYNThiApp8jTgeHzDzIbB5DCI= 3675109    0s         9209444

$ # verify
$ echo 4JFWFevraBg5yLD2rCFzGMpbKZOC3BV5rNom+Um7EGg= | wg pubkey 
2025Q6KcDb+v/nj2/ErYNThiApp8jTgeHzDzIbB5DCI=
```

or run the tool from the source repository:
```console
$ go run . --prefix=2025
```

or use Docker image:
```console
$ docker pull ghcr.io/alexanderyastrebov/wireguard-vanity-key:latest
$ docker run ghcr.io/alexanderyastrebov/wireguard-vanity-key:latest --prefix=2025
```

## Benchmark

The tool checks ~12'000'000 keys per second on a test machine:

```console
$ go test . -run=NONE -bench=BenchmarkFindPointParallel -benchmem -count=10
goos: linux
goarch: amd64
pkg: github.com/AlexanderYastrebov/wireguard-vanity-key
cpu: Intel(R) Core(TM) i5-8350U CPU @ 1.70GHz
BenchmarkFindPointParallel-8    13230948                86.06 ns/op            0 B/op          0 allocs/op
BenchmarkFindPointParallel-8    13185460                87.24 ns/op            0 B/op          0 allocs/op
BenchmarkFindPointParallel-8    12931213                88.63 ns/op            0 B/op          0 allocs/op
BenchmarkFindPointParallel-8    12287206                89.30 ns/op            0 B/op          0 allocs/op
BenchmarkFindPointParallel-8    12603378                90.61 ns/op            0 B/op          0 allocs/op
BenchmarkFindPointParallel-8    12401572                91.26 ns/op            0 B/op          0 allocs/op
BenchmarkFindPointParallel-8    12177254                93.40 ns/op            0 B/op          0 allocs/op
BenchmarkFindPointParallel-8    12052425                92.88 ns/op            0 B/op          0 allocs/op
BenchmarkFindPointParallel-8    12028226                93.76 ns/op            0 B/op          0 allocs/op
BenchmarkFindPointParallel-8    12186556                93.87 ns/op            0 B/op          0 allocs/op
PASS
ok      github.com/AlexanderYastrebov/wireguard-vanity-key      21.154s
```

## Similar projects

* [wireguard-vanity-address](https://github.com/warner/wireguard-vanity-address)
* [wireguard-vanity-keygen](https://github.com/axllent/wireguard-vanity-keygen)
* [Wireguard-Vanity-Key-Searcher](https://github.com/volleybus/Wireguard-Vanity-Key-Searcher)
* [wgmine](https://github.com/thatsed/wgmine)
* [Vanity](https://github.com/samuel-lucas6/Vanity)

# wireguard-vanity-key

Inspired by [wireguard-vanity-address "faster algorithm"](https://github.com/warner/wireguard-vanity-address/pull/15),
this tool searches for a [WireGuard](https://www.wireguard.com/) Curve25519 keypair
with a base64-encoded public key that has a specified prefix.

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

or use Docker image:
```console
$ docker pull ghcr.io/alexanderyastrebov/wireguard-vanity-key:latest
$ docker run  ghcr.io/alexanderyastrebov/wireguard-vanity-key:latest -prefix=2025
```

## Performance

The tool checks ~18'000'000 keys per second on a test machine:

```console
$ go test . -run=NONE -bench=BenchmarkFindPointParallel -benchmem -count=10
goos: linux
goarch: amd64
pkg: github.com/AlexanderYastrebov/wireguard-vanity-key
cpu: Intel(R) Core(TM) i5-8350U CPU @ 1.70GHz
BenchmarkFindPointParallel-8    19739348                54.33 ns/op            0 B/op          0 allocs/op
BenchmarkFindPointParallel-8    19185619                55.42 ns/op            0 B/op          0 allocs/op
BenchmarkFindPointParallel-8    19316592                56.08 ns/op            0 B/op          0 allocs/op
BenchmarkFindPointParallel-8    18855543                56.95 ns/op            0 B/op          0 allocs/op
BenchmarkFindPointParallel-8    18705961                56.46 ns/op            0 B/op          0 allocs/op
BenchmarkFindPointParallel-8    18718236                56.45 ns/op            0 B/op          0 allocs/op
BenchmarkFindPointParallel-8    18693268                56.78 ns/op            0 B/op          0 allocs/op
BenchmarkFindPointParallel-8    18495776                57.85 ns/op            0 B/op          0 allocs/op
BenchmarkFindPointParallel-8    18160232                58.81 ns/op            0 B/op          0 allocs/op
BenchmarkFindPointParallel-8    18100197                57.53 ns/op            0 B/op          0 allocs/op
PASS
ok      github.com/AlexanderYastrebov/wireguard-vanity-key      21.154s
```

Time to get n-symbol prefix at rate 18'000'000 keys per second with probability of:
```
n                  50%                  95%                  99%
4              0:00:01              0:00:03              0:00:04
5              0:00:41              0:02:59              0:04:35
6              0:44:06              3:10:37              4:53:01
7      1 day, 23:02:41     8 days, 11:19:25     13 days, 0:33:28
8   125 days, 10:51:28    542 days, 4:42:40   833 days, 11:42:23
```
I.e. you may expect a match for a 5 character prefix after
41 seconds 50% of the time and after 4-5 minutes to be 99% sure.

## Similar projects

* [wireguard-vanity-address](https://github.com/warner/wireguard-vanity-address)
* [wireguard-vanity-keygen](https://github.com/axllent/wireguard-vanity-keygen)
* [Wireguard-Vanity-Key-Searcher](https://github.com/volleybus/Wireguard-Vanity-Key-Searcher)
* [wgmine](https://github.com/thatsed/wgmine)
* [Vanity](https://github.com/samuel-lucas6/Vanity)
* [mkp224o](https://github.com/cathugger/mkp224o)

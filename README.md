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

The tool checks ~29'000'000 keys per second on a test machine:

```console
$ go run . -prefix=GoodLuckWithThisPrefix -timeout=20s
private                                      public                                       attempts   duration   attempts/s
-                                            GoodLuckWithThisPrefix...                    583920640  20s        29194831
```

In practice it finds 4-character prefix in a second and 5-character prefix in a minute:
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

The tool supports blind search, i.e. when worker does not know the private key, see [demo-blind.sh](demo-blind.sh).

## Similar projects

* [wireguard-vanity-address](https://github.com/warner/wireguard-vanity-address)
* [wireguard-vanity-keygen](https://github.com/axllent/wireguard-vanity-keygen)
* [Wireguard-Vanity-Key-Searcher](https://github.com/volleybus/Wireguard-Vanity-Key-Searcher)
* [wgmine](https://github.com/thatsed/wgmine)
* [Vanity](https://github.com/samuel-lucas6/Vanity)

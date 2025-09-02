# wireguard-vanity-key

This tool searches for a [WireGuard](https://www.wireguard.com/) Curve25519 keypair
with a base64-encoded public key that has a specified prefix.

Compared to [similar tools](#similar-tools), it uses [the fastest search algorithm](#the-fastest-search-algorithm) 🚀

## Example

Install the tool locally and run:
```console
$ go install github.com/AlexanderYastrebov/wireguard-vanity-key@latest
$ wireguard-vanity-key --prefix=2025
private                                      public                                       attempts   duration   attempts/s
MuWXIBrtHVtgbTUe+1ujcaJCN6P/NnISk25dD85XsBs= 2025eq2dnBEQ3A7ak0eSwND4U9JxjNKUQ4m7407P2hE= 14478798   0s         41972460

$ # verify
$ echo MuWXIBrtHVtgbTUe+1ujcaJCN6P/NnISk25dD85XsBs= | wg pubkey
2025eq2dnBEQ3A7ak0eSwND4U9JxjNKUQ4m7407P2hE=
```

or run the tool from the source repository:
```console
$ go run . --prefix=2025
```

or use the Docker image:
```console
$ docker pull ghcr.io/alexanderyastrebov/wireguard-vanity-key:latest
$ docker run  ghcr.io/alexanderyastrebov/wireguard-vanity-key:latest --prefix=2025
```

## Performance

The tool checks ~65'000'000 keys per second on a test machine:

```console
$ go run . --prefix=GoodLuckWithThisPrefix --timeout=20s
private                                      public                                       attempts   duration   attempts/s
-                                            GoodLuckWithThisPrefix...                    1374379620 20s        68701706
```

In practice, it finds a 4-character prefix in a second and a 5-character prefix in a minute:
```console
$ while go run . --prefix=AYAYA ; do : ; done
private                                      public                                       attempts   duration   attempts/s
P15GaB7DdiKzYfCZynsClGhHVuIFilvV3zaNCgRuPn4= AYAYA++HaWPnqaCBE+xlB+GNvyUI3LLtLAoHhMIEBgs= 219066590  3s         69609692
private                                      public                                       attempts   duration   attempts/s
7pqZA+StBeqiLjltLJQhzvCzuj75tx6HiRAOxOF68GE= AYAYAGeWrXfsnvbsv1n7FkZG3zgGmYc6sjr2HUxLJXI= 530717186  8s         67577159
private                                      public                                       attempts   duration   attempts/s
M2bhXvumSXmjiiHxgV7u8uBfbS9hTXdnXlrZUYNZ1gc= AYAYAPY8KNO/tm6vLITC32eT/IKn2V/lVIz+0lcab1I= 3458371931 58s        59901869
private                                      public                                       attempts   duration   attempts/s
ALnF+zYUI9ZCGCP3/PwgkBQrVNHIsXJdZ+AYEbYt9BU= AYAYAZd81yqulNeXXi23mHIMQp8q/MecGNxUOkhXdiI= 1197282698 21s        57711372
private                                      public                                       attempts   duration   attempts/s
koJOikpom8J7fEw1E7Gb3TU5A0PBGtqGXP9KA09Mtm8= AYAYANW4FVp515jAwSpZcphIzbfjKAeuzvZVTyiERSk= 133529424  2s         55282205
```

Each additional character increases search time by a factor of 64.

## Blind search

The tool supports blind search, i.e., when the worker does not know the private key. See [demo-blind.sh](demo-blind.sh).

## Kubernetes

You can run the tool in a distributed manner in Kubernetes cluster using the [demo-k8s.yaml](demo-k8s.yaml) manifest
to search for a vanity key without exposing the private key:

```console
$ # Generate secure starting key pair
$ wg genkey | tee /dev/stderr | wg pubkey
YI5+UcKmyLdeRDqU8l3k53wrUZO9Mw23NpvB8tDtvWU=
startkQgqI9Gv1IX7eNa2qeFhpYBRDwpz40JIAAYOSk=

$ # Edit demo-k8s.yaml to configure prefix, starting public key, parallelism, and resource limits 💸

$ # Create search job
$ kubectl apply -f demo-k8s.yaml
job.batch/wvk created

$ # Check job
$ kubectl get job wvk
NAME   STATUS    COMPLETIONS   DURATION   AGE
wvk    Running   0/10          2m53s      2m53s

$ # Check pods
$ kubectl get pods --selector=batch.kubernetes.io/job-name=wvk
NAME         READY   STATUS    RESTARTS   AGE
wvk-0-8tdz5  1/1     Running   0          3m8s
wvk-1-pmnkn  1/1     Running   0          3m8s
wvk-2-2ls7m  1/1     Running   0          3m8s
wvk-3-rd7gx  1/1     Running   0          3m8s
wvk-4-jqksz  1/1     Running   0          3m8s
wvk-5-vj6gd  1/1     Running   0          3m8s
wvk-6-vhgmc  1/1     Running   0          3m8s
wvk-7-drr98  1/1     Running   0          3m8s
wvk-8-tmb6c  1/1     Running   0          3m8s
wvk-9-gxlp2  1/1     Running   0          3m8s

$ # Check resource usage
$ kubectl top pods --selector=batch.kubernetes.io/job-name=wvk

$ # Wait for the job to complete
$ kubectl wait --for=condition=complete job/wvk --timeout=1h
job.batch/wvk condition met

$ # Job is complete
$ kubectl get job wvk
NAME   STATUS     COMPLETIONS   DURATION   AGE
wvk    Complete   1/999999      34m        37m

$ # Get found offset from the logs
$ kubectl logs jobs/wvk
7538451707115552752

$ # Generate new private vanity key by offsetting the starting private key
$ echo YI5+UcKmyLdeRDqU8l3k53wrUZO9Mw23NpvB8tDtvWU= | wireguard-vanity-key add --offset=7538451707115552752
4I4EWan32HJbRDqU8l3k53wrUZO9Mw23NpvB8tDtvWU=

$ # Get the vanity public key
$ echo 4I4EWan32HJbRDqU8l3k53wrUZO9Mw23NpvB8tDtvWU= | wg pubkey
wvk+k8shgsJcW5EKet2AkViKc7a/0Ud8/EDOy91aCQg=

$ # Delete the job
$ kubectl delete job wvk
job.batch "wvk" deleted
```

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

This tool uses only **3.5 (three and a half)** field multiplications per candidate key and other optimizations, which makes it the fastest 🚀

See [vanity25519](https://github.com/AlexanderYastrebov/vanity25519) for algorithm implementation.

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

### 🧮 Batch field inversion

Getting the public key in WireGuard format (Montgomery form) requires a field inversion, which uses **265** field multiplications.
This tool uses the Montgomery trick to implement batch inversion, so for a batch of N candidates,
only **1 (one)** inversion is needed, resulting in a huge speedup.

### 🪞 Montgomery coordinates and offset symmetry

The algorithm uses Montgomery coordinates and exploits symmetries in precomputed point offsets,
saving even more field multiplications.

### ⚡ Fast base64 prefix check

Other tools encode the full public key to base64 and compare the prefix. This tool decodes the base64 prefix and compares it to public key bytes directly. See https://github.com/AlexanderYastrebov/wireguard-vanity-key/pull/5.

### 🏆 High-performance C implementation

For raw speed, the C worker (`wvk`) uses [awslabs/s2n-bignum](https://github.com/awslabs/s2n-bignum) -
a highly optimized field arithmetic library written in assembly.
The worker supports prefix lengths up to 10 base64 characters, so the prefix check becomes a single masked integer comparison.
These two optimizations make `wvk` ~2 times faster than the Go implementation.
See https://github.com/AlexanderYastrebov/wireguard-vanity-key/pull/15.

# Set ARCH to "x86" or "arm" as appropriate; this is literally the
# directory name in which the libs2nbignum.a file should be found.
# If this block fails to do it automatically, replace with manual
# setting: ARCH=x86 or ARCH=arm as appropriate

ARCH=unknown_architecture

UNAME_RESULT=$(shell uname -m)

ifeq ($(UNAME_RESULT),x86_64)
ARCH=x86
endif

ifeq ($(UNAME_RESULT),aarch64)
ARCH=arm
endif

ifeq ($(UNAME_RESULT),arm64)
ARCH=arm
endif

wvk: s2nbignum main.c
	# Use -pg to profile
	gcc -static -Wall -O3 -o wvk main.c -I./s2n-bignum/include/ -L./s2n-bignum/$(ARCH)/ -ls2nbignum

.PHONY: test
test: s2nbignum main_test.c
	gcc -o wvk_test -nostartfiles -Wl,--entry=test_all main_test.c -I./s2n-bignum/include/ -L./s2n-bignum/x86/ -ls2nbignum
	./wvk_test

main.S: s2nbignum main.c
	gcc -O3 -S -o main.S main.c -I./s2n-bignum/include/ -L./s2n-bignum/$(ARCH)/ -ls2nbignum

.PHONY: s2nbignum
s2nbignum:
	if [ ! -d ./s2n-bignum/ ]; then \
		git clone --depth=1 git@github.com:awslabs/s2n-bignum.git; \
	fi
	make -C ./s2n-bignum/$(ARCH)/

.PHONY: clean
clean:
	rm -f wvk wvk_test main.S
	make -C ./s2n-bignum/$(ARCH)/ clean

.PHONY: benchmark
benchmark: wvk
	timeout -s SIGINT 20s ./wvk offset AYAYAiziAiJvySpk0uHM6qNZgUYLDNxqYNnyI3v9mRQ= AYAYAYAYA 0 0 || exit 0

.PHONY: benchmark-parallel
benchmark-parallel: wvk
	timeout -s SIGINT 20s ./pwvk.sh AYAYAYAYA || exit 0

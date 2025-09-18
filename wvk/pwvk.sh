#!/bin/bash
#
# This script searches for public key prefix using multiple worker processes.
#
set -euo pipefail

# The prefix to search for
prefix=${1:-"AY/"}

# Generate secure staring private key
private=$(wg genkey)

# and its public key
public=$(echo $private | wg pubkey)

# Launch parallel workers and capture their stdout and stderr
wstdout=$(mktemp)
wstderr=$(mktemp)

for i in $(seq 1 $(nproc))
do
    echo $((RANDOM*$RANDOM*$RANDOM*$RANDOM))
done |
parallel --ungroup --halt now,success=1 --termseq SIGINT,100 \
    ./wvk offset $public $prefix {} 0 \
    1>$wstdout 2>$wstderr &

ppid=$!

# Print sum of attempts/s
function print_attempts()
{
    attempts=$(( $(cat $wstderr | grep -F 'attempts/s' | sed 's#attempts/s: ##' | tr "\n" "+" ) 0 ))
    echo "attempts/s: $attempts" >&2
}

# Handle signals and exit
function interrupt()
{
    kill -s SIGTERM $ppid
    wait
    print_attempts
    exit 1
}
trap interrupt SIGINT SIGTERM

wait

print_attempts

# Get the first found offset
offset=$(head -n 1 $wstdout)

# Generate new private vanity key by offsetting the starting $private key.
private_vanity=$(echo $private | ./wvk add $offset)

# Print vanity key pair
echo $private_vanity
echo $private_vanity | wg pubkey

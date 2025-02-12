#!/bin/bash
#
# This script demonstrates blind prefix search,
# i.e. when worker does not know the private key.
#

prefix=AY/

# Generate secure staring private key
private=$(wg genkey)

# and its public key
public=$(echo $private | wg pubkey)

# Search for prefix by incrementing $public key and report result offset.
# Note that $private key is not involved and search can be scaled horizontally.
offset=$(wireguard-vanity-key --prefix=$prefix --public=$public --output=offset)

# Generate new private vanity key by offsetting the starting $private key.
# Give it the $prefix to choose from two possible results.
private_vanity=$(echo $private | wireguard-vanity-key add --offset=$offset --prefix=$prefix)

# Print vanity key pair
echo $private_vanity
echo $private_vanity | wg pubkey

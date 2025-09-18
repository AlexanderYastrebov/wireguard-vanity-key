#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <signal.h>
#include <inttypes.h>
#include <byteswap.h>
#include <time.h>
#include "s2n-bignum.h"

#define EXIT_USAGE 2
#define EXIT_INTERRUPTED 3

int usage()
{
    fprintf(stderr, "Usage: wvk offset PUBLIC_KEY PREFIX SKIP LIMIT\n");
    fprintf(stderr, "Find PREFIX offset from PUBLIC_KEY after skipping SKIP steps and making up to LIMIT steps\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Usage: wvk add OFFSET\n");
    fprintf(stderr, "Add OFFSET to the private key read from stdin\n");
    return EXIT_USAGE;
}

static const int batch_size = 4 * 1024; // must be even

typedef uint64_t field_element[4];

static const int X = 0;
static const int Y = 1;
static const int Z = 2;
static const int T = 3;

typedef uint64_t edwards_point[8];
typedef field_element edwards_point_ep[4]; // X, Y, Z, T
typedef field_element montgomery_point[2]; // X, Y

void make_offsets(montgomery_point *offsets, int n);
void invert_batch(field_element *a, field_element *b, int n);
int montgomery_from_bytes_with_offset(montgomery_point result, const uint8_t *public_key, uint64_t offset);
void montgomery_add_batch_x(field_element *x, montgomery_point p1, montgomery_point *p2s, field_element *dx, int n);

void montgomery_from_edwards(montgomery_point m, edwards_point p);
void montgomery_from_edwards_ep(montgomery_point m, edwards_point_ep p);
void montgomery_set(montgomery_point dst, montgomery_point src);
void montgomery_add(montgomery_point result, montgomery_point p1, montgomery_point p2);
void montgomery_add_dxinv(montgomery_point result, montgomery_point p1, montgomery_point p2, field_element *dxinv);

void edwards_bytes_montgomery(edwards_point p, uint8_t *dst);
void edwards_from_edwards_ep(edwards_point v, edwards_point_ep p);
void edwards_ep_from_edwards(edwards_point_ep v, edwards_point p);

void field_element_set(field_element v, field_element x);
void field_element_print_base64(field_element x);

uint64_t reverse_bits(uint64_t x);

int base64_encode(const uint8_t *src, int srclen, char *dst);
int base64_decode(const char *src, int srclen, uint8_t *dst);

static field_element ONE = {
    UINT64_C(1),
    UINT64_C(0),
    UINT64_C(0),
    UINT64_C(0),
};

static field_element COFACTOR = {
    UINT64_C(8),
    UINT64_C(0),
    UINT64_C(0),
    UINT64_C(0),
};

// Constant -|sqrt(-486664)| - the scaling factor for bi-rational map
// calculated such that Edwards generator point maps to Montgomery base point.
static field_element SQRT486664 = {
    UINT64_C(3716027510060384743),
    UINT64_C(4205847681119217021),
    UINT64_C(3280018162556579969),
    UINT64_C(8131550443321948484),
};

// Montgomery "curve25519" v^2 = u^3 + A*u^2 + u parameter
// Constant A = 486662
static field_element A = {
    UINT64_C(486662),
    UINT64_C(0),
    UINT64_C(0),
    UINT64_C(0),
};

static uint8_t prefix_bytes[32];
static int prefix_bits;

sig_atomic_t interrupted = 0;

void signal_handler(int sig)
{
    interrupted = 1;
}

int cmd_add(int argc, char *argv[]);

int main(int argc, char *argv[])
{
    if (argc > 1 && !strcmp("add", argv[1]))
    {
        return cmd_add(argc - 1, argv + 1);
    }
    if (argc != 6 || strcmp("offset", argv[1]))
    {
        return usage();
    }
    const char *arg_public_key = argv[2];
    const char *arg_prefix = argv[3];
    const char *arg_skip = argv[4];
    const char *arg_limit = argv[5];

    signal(SIGINT, signal_handler);

    uint8_t public_key[32];
    {
        if (strlen(arg_public_key) == 44)
        {
            base64_decode(arg_public_key, 44, public_key);
        }
        else
        {
            fprintf(stderr, "Invalid public key\n");
            return usage();
        }
    }

    uint64_t skip = strtoul(arg_skip, NULL, 10);
    uint64_t limit = strtoul(arg_limit, NULL, 10);

    uint64_t prefix_match, mask;
    {
        int prefix_len = strlen(arg_prefix);
        prefix_bits = 6 * prefix_len;

        // limit prefix length to 64 bits for fast testing via uint64_t mask
        if (prefix_bits > 64)
        {
            fprintf(stderr, "Maximum supported prefix length is 64 bits ~ 10 base64 characters\n");
            return usage();
        }
        base64_decode(arg_prefix, prefix_len, prefix_bytes);

        mask = reverse_bits(bswap_64((1ul << prefix_bits) - 1));

        uint64_t t[4];
        bignum_fromlebytes_4(t, prefix_bytes);
        prefix_match = t[0] & mask;
    }

    montgomery_point p;
    {
        if (montgomery_from_bytes_with_offset(p, public_key, skip) != 0)
        {
            fprintf(stderr, "Invalid public key\n");
            return EXIT_FAILURE;
        }
    }

    field_element x[batch_size];
    field_element dx[batch_size / 2 + 1];
    montgomery_point offsets[batch_size / 2];
    montgomery_point batch_offset;

    make_offsets(offsets, batch_size / 2);

    montgomery_set(batch_offset, offsets[0]);
    montgomery_add(batch_offset, batch_offset, offsets[batch_size / 2 - 1]);
    montgomery_add(batch_offset, batch_offset, offsets[batch_size / 2 - 1]);

    // Shift start point by half batch to avoid negative offsets
    montgomery_add(p, p, offsets[batch_size / 2 - 1]);
    uint64_t n = batch_size / 2;

    clock_t start_time = clock();
    while (!interrupted)
    {
        // montgomery_add_batch_x inverts last element of dx
        bignum_sub_p25519(dx[batch_size / 2], batch_offset[X], p[X]);

        montgomery_add_batch_x(x, p, offsets, dx, batch_size / 2);

        for (int i = 0; i < batch_size; i++)
        {
            if ((x[i][0] & mask) == prefix_match)
            {
                if (i < batch_size / 2)
                {
                    n += (i + 1);
                }
                else
                {
                    n -= (i + 1 - batch_size / 2);
                }
                goto found;
            }
        }

        if ((p[X][0] & mask) == prefix_match)
        {
            goto found;
        }

        // Advance to next batch: p = p + batch_offset
        // dx[batch_size/2] was calculated above and inverted in montgomery_add_batch_x
        montgomery_add_dxinv(p, p, batch_offset, &dx[batch_size / 2]);

        n += batch_size + 1;

        if (limit > 0)
        {
            if (limit <= batch_size + 1)
            {
                break;
            }
            limit -= batch_size + 1;
        }
    }
    double seconds;
found:
    seconds = (double)(clock() - start_time) / (double)CLOCKS_PER_SEC;

    printf("%lu\n", skip + n);

    fprintf(stderr, "seconds: %0.0f\n", seconds);
    fprintf(stderr, "attempts/s: %0.0f\n", n / seconds);

    return interrupted ? EXIT_INTERRUPTED : EXIT_SUCCESS;
}

// add OFFSET
int cmd_add(int argc, char *argv[])
{
    const char *arg_offset = argv[1];

    field_element s0;
    {
        uint8_t buf32[32];
        char buf[44];
        if (fgets(buf, 44, stdin) == NULL)
        {
            fprintf(stderr, "Failed to read private key from stdin\n");
            usage();
        }
        base64_decode(buf, 44, buf32);
        bignum_fromlebytes_4(s0, buf32);
    }

    uint64_t offset = strtoul(arg_offset, NULL, 10);

    // Get start public key from start private key
    uint8_t start_public_key[32];
    {
        edwards_point p;
        edwards25519_scalarmulbase(p, s0);
        edwards_bytes_montgomery(p, start_public_key);
    }

    // Calculate expected vanity public key by adding offset to start public key
    uint8_t vanity_public_key[32];
    {
        montgomery_point p;
        if (montgomery_from_bytes_with_offset(p, start_public_key, offset) != 0)
        {
            fprintf(stderr, "Invalid start public key\n");
            return EXIT_FAILURE;
        }
        bignum_tolebytes_4(vanity_public_key, p[X]);
    }

    field_element so = {offset, 0, 0, 0};
    bignum_mul_p25519(so, so, COFACTOR);

    // Try both s0 + so and s0 - so
    field_element sp, sm;
    bignum_add_p25519(sp, s0, so);
    bignum_sub_p25519(sm, s0, so);

    edwards_point pp, pm;
    edwards25519_scalarmulbase(pp, sp);
    edwards25519_scalarmulbase(pm, sm);

    uint8_t buf32[32];

    edwards_bytes_montgomery(pp, buf32);
    if (!memcmp(buf32, vanity_public_key, 32))
    {
        field_element_print_base64(sp);
        return EXIT_SUCCESS;
    }

    edwards_bytes_montgomery(pm, buf32);
    if (!memcmp(buf32, vanity_public_key, 32))
    {
        field_element_print_base64(sm);
        return EXIT_SUCCESS;
    }

    fprintf(stderr, "Offset does not match private key\n");
    return EXIT_FAILURE;
}

// Generate Montgomery offset points: offsets[i] = B*8*(i+1).
void make_offsets(montgomery_point *offsets, int n)
{
    edwards_point p;
    for (int i = 0; i < n; i++)
    {
        field_element oi = {(i + 1), 0, 0, 0};
        bignum_mul_p25519(oi, oi, COFACTOR);
        edwards25519_scalarmulbase(p, oi);
        montgomery_from_edwards(offsets[i], p);
    }
}

// invert_batch calculates a[i] = 1/a[i] using b as a scratch buffer.
//
// It uses:
//
//	3*(n-1) multiplications
//	1 invert = ~265 multiplications
//
// Complexity: 3M*n + 262M
//
// https://en.wikipedia.org/wiki/Modular_multiplicative_inverse#Multiple_inverses
void invert_batch(field_element *a, field_element *b, int n)
{
    field_element t;
    field_element pa; // product a[0]*a[1]*...*a[i]

    field_element_set(pa, a[0]);
    for (int i = 1; i < n; i++)
    {
        field_element_set(b[i], pa);
        bignum_mul_p25519(pa, pa, a[i]);
    }

    field_element paInv;
    bignum_inv_p25519(paInv, pa);

    for (int i = n - 1; i > 0; i--)
    {
        bignum_mul_p25519(t, paInv, b[i]);
        bignum_mul_p25519(paInv, paInv, a[i]);
        field_element_set(a[i], t);
    }
    field_element_set(a[0], paInv);
}

// Converts Montgomery u-coordinate bytes to Montgomery point and adds B*8*offset
// Returns 0 on success, -1 on error
int montgomery_from_bytes_with_offset(montgomery_point result, const uint8_t *public_key, uint64_t offset)
{
    field_element u, y, t;
    edwards_point tp;
    edwards_point_ep p, po;
    uint8_t buf[32];

    // Convert bytes to field element
    bignum_fromlebytes_4(u, public_key);

    // Convert Montgomery u-coordinate to Edwards y-coordinate
    // y = (u - 1) / (u + 1)
    bignum_add_p25519(t, u, ONE);
    bignum_inv_p25519(t, t);
    bignum_sub_p25519(y, u, ONE);
    bignum_mul_p25519(y, y, t);

    // Convert to bytes and decode as Edwards point
    bignum_tolebytes_4(buf, y);
    if (edwards25519_decode(tp, buf) != 0)
    {
        return -1; // Invalid point
    }
    edwards_ep_from_edwards(p, tp);

    field_element so = {offset, 0, 0, 0};
    bignum_mul_p25519(so, so, COFACTOR);
    edwards25519_scalarmulbase(tp, so);
    edwards_ep_from_edwards(po, tp);

    edwards25519_epadd((uint64_t *)p, (uint64_t *)p, (uint64_t *)po);

    montgomery_from_edwards_ep(result, p);

    return 0;
}

// montgomery_add_batch_x adds a batch of n points to a given point p1 and
// returns x-coordinates of the 2*n resulting points:
//
//	p3x = {p1 + p2s[0], ... , p1 + p2s[n-1], p1 - p2s[0], ... , p1 - p2s[n-1]}
//
// Complexity for 2*n resulting x-coordinates:
//
//	(4M + 6A)*n + 3M*(n+1) + 262M + 1A = (7M + 6A)*n + 265M + 1A
//
// I.e. (3.5M + 3A) per resulting x-coordinate amortized.
//
// It requires:
//
//	len(dx) = n+1, uses dx[:n] as a scratch buffer and calculates dx[n] = 1/dx[n]
//	len(p3x)  = 2*n, for resulting x-coordinates
void montgomery_add_batch_x(field_element *p3x, montgomery_point p1, montgomery_point *p2s, field_element *dx, int n)
{
    field_element Ax1, Ax1x2, t;

    for (int i = 0; i < n; i++)
    {
        bignum_sub_p25519(dx[i], p2s[i][X], p1[X]);
    }

    // dx has extra last element to invert
    // use p3x as a scratch buffer
    invert_batch(dx, p3x, n + 1);

    bignum_add_p25519(Ax1, A, p1[X]);

    for (int i = 0; i < n; i++)
    {
        field_element *p2 = p2s[i];

        bignum_add_p25519(Ax1x2, Ax1, p2[X]);

        // For p1 + p2: slope = (y2 - y1) / (x2 - x1)
        bignum_sub_p25519(t, p2[Y], p1[Y]);
        bignum_mul_p25519(t, t, dx[i]);
        bignum_sqr_p25519(t, t);
        bignum_sub_p25519(p3x[i], t, Ax1x2);

        // For p1 - p2: slope = (-y2 - y1) / (x2 - x1) = (y1 + y2) / (x1 - x2)
        bignum_add_p25519(t, p2[Y], p1[Y]);
        bignum_mul_p25519(t, t, dx[i]);
        bignum_sqr_p25519(t, t);
        bignum_sub_p25519(p3x[n + i], t, Ax1x2);
    }
}

// Montgomery point arithmetic functions

// https://www.rfc-editor.org/rfc/rfc7748.html#section-4.1
// (u, v) = ((1+y)/(1-y), sqrt(-486664)*u/x)
void montgomery_from_edwards(montgomery_point m, edwards_point p)
{
    field_element t;
    uint64_t *x = p;
    uint64_t *y = p + 4;
    uint64_t *u = m[X];
    uint64_t *v = m[Y];

    bignum_sub_p25519(t, ONE, y);
    bignum_inv_p25519(t, t);
    bignum_add_p25519(u, ONE, y);
    bignum_mul_p25519(u, u, t); // u = (1+y)/(1-y)

    bignum_inv_p25519(t, x);
    bignum_mul_p25519(v, SQRT486664, u);
    bignum_mul_p25519(v, v, t); // v = sqrt(-486664)*u/x
}

void montgomery_from_edwards_ep(montgomery_point m, edwards_point_ep p)
{
    edwards_point pe;
    edwards_from_edwards_ep(pe, p);
    montgomery_from_edwards(m, pe);
}

void montgomery_set(montgomery_point dst, montgomery_point src)
{
    field_element_set(dst[X], src[X]);
    field_element_set(dst[Y], src[Y]);
}

// https://en.wikipedia.org/wiki/Montgomery_curve
//
// Montgomery curve point addition formulae:
// x3 = ((y2 - y1) / (x2 - x1))^2 - A - x1 - x2
// y3 = (2*x1 + x2 + A) * ((y2 - y1) / (x2 - x1)) - ((y2 - y1) / (x2 - x1))^3 - y1
//
// Complexity: 1I + 4M + 7A
void montgomery_add(montgomery_point result, montgomery_point p1, montgomery_point p2)
{
    field_element dxInv;
    bignum_sub_p25519(dxInv, p2[X], p1[X]);
    bignum_inv_p25519(dxInv, dxInv);

    montgomery_add_dxinv(result, p1, p2, &dxInv);
}

// Montgomery curve point addition with precomputed dxinv = 1/(x2 - x1)
// Complexity: 4M + 7A
void montgomery_add_dxinv(montgomery_point p3, montgomery_point p1, montgomery_point p2, field_element *dxinv)
{
    field_element x3, y3;
    field_element dy, slope, slopeSquared, slopeCubed, x2A;
    field_element x12A, xSum;

    // x2A = x2 + A
    bignum_add_p25519(x2A, p2[X], A);

    // dy = y2 - y1
    bignum_sub_p25519(dy, p2[Y], p1[Y]);

    // slope = dy * dxinv
    bignum_mul_p25519(slope, dy, *dxinv);
    bignum_mul_p25519(slopeSquared, slope, slope);
    bignum_mul_p25519(slopeCubed, slopeSquared, slope);

    // x12A = x1 + x2 + A
    bignum_add_p25519(x12A, p1[X], x2A);

    // x3 = slope^2 - x1 - x2 - A
    bignum_sub_p25519(x3, slopeSquared, x12A);

    // xSum = 2*x1 + x2 + A
    bignum_add_p25519(xSum, p1[X], x12A);

    // y3 = xSum * slope - slope^3 - y1
    bignum_mul_p25519(y3, xSum, slope);
    bignum_sub_p25519(y3, y3, slopeCubed);
    bignum_sub_p25519(y3, y3, p1[Y]);

    field_element_set(p3[X], x3);
    field_element_set(p3[Y], y3);
}

void edwards_bytes_montgomery(edwards_point p, uint8_t *dst)
{
    // u = (1 + y) / (1 - y)
    uint64_t *y = p + 4;
    field_element n, d, t, u;

    bignum_add_p25519(n, ONE, y);
    bignum_sub_p25519(d, ONE, y);
    bignum_inv_p25519(t, d);
    bignum_mul_p25519(u, n, t);

    bignum_tolebytes_4(dst, u);
}

void edwards_from_edwards_ep(edwards_point v, edwards_point_ep p)
{
    field_element t;
    bignum_inv_p25519(t, p[Z]);
    bignum_mul_p25519(v, p[X], t);     // x = X/Z
    bignum_mul_p25519(v + 4, p[Y], t); // y = Y/Z
}

void edwards_ep_from_edwards(edwards_point_ep v, edwards_point p)
{
    field_element_set(v[X], p);
    field_element_set(v[Y], p + 4);
    field_element_set(v[Z], ONE);
    bignum_mul_p25519(v[T], v[X], v[Y]);
}

void field_element_set(field_element v, field_element x)
{
    v[0] = x[0];
    v[1] = x[1];
    v[2] = x[2];
    v[3] = x[3];
}

void field_element_print_base64(field_element x)
{
    uint8_t buf32[32];
    char buf45[45];
    buf45[44] = 0;

    bignum_tolebytes_4(buf32, x);
    base64_encode(buf32, 32, buf45);
    printf("%s\n", buf45);
}

uint64_t reverse_bits(uint64_t x)
{
    uint64_t result = 0;
    for (int i = 0; i < 64; i++)
    {
        result <<= 1;
        result |= x & 1;
        x >>= 1;
    }
    return result;
}

static const char *base64_digits = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";

int base64_decode(const char *src, int srclen, uint8_t *dst)
{
    int j = 0;
    uint32_t s = 0;
    int bits = 0;

    for (int i = 0; i < srclen && src[i] != '='; i++)
    {
        char *p = strchr(base64_digits, src[i]);
        if (p == NULL)
        {
            return -1;
        }
        int d = p - base64_digits;

        s <<= 6;
        s |= d;
        bits += 6;

        if (bits == 24)
        {
            dst[j++] = (s >> 16) & 0xff;
            dst[j++] = (s >> 8) & 0xff;
            dst[j++] = (s >> 0) & 0xff;

            bits = 0;
            s = 0;
        }
    }

    switch (bits)
    {
    case 0:
        break;
    case 18:
        dst[j++] = (s >> 10) & 0xff;
        dst[j++] = (s >> 2) & 0xff;
        dst[j++] = (s << 6) & 0xff;
        break;
    case 12:
        dst[j++] = (s >> 4) & 0xff;
        dst[j++] = (s << 4) & 0xff;
        break;
    case 6:
        dst[j++] = (s << 2) & 0xff;
        break;
    default:
        return -1;
    }
    return j;
}

int base64_encode(const uint8_t *src, int srclen, char *dst)
{
    int j = 0;
    int i = 0;
    uint32_t s = 0;

    for (; i < (srclen / 3) * 3; i += 3)
    {
        s = src[i] << 16;
        s |= src[i + 1] << 8;
        s |= src[i + 2];

        dst[j++] = base64_digits[(s >> 18) & 0b111111];
        dst[j++] = base64_digits[(s >> 12) & 0b111111];
        dst[j++] = base64_digits[(s >> 6) & 0b111111];
        dst[j++] = base64_digits[(s >> 0) & 0b111111];
    }

    switch (srclen - i)
    {
    case 0:
        break;
    case 1:
        dst[j++] = base64_digits[(src[i] >> 2) & 0b111111];
        dst[j++] = base64_digits[(src[i] << 4) & 0b111111];
        dst[j++] = '=';
        dst[j++] = '=';
        break;
    case 2:
        dst[j++] = base64_digits[(src[i] >> 2) & 0b111111];
        dst[j++] = base64_digits[((src[i] << 4) & 0b110000) | ((src[i + 1] >> 4) & 0b001111)];
        i++;
        dst[j++] = base64_digits[(src[i] << 2) & 0b111111];
        dst[j++] = '=';
        break;
    default:
        return -1;
    }
    return j;
}
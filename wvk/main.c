#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <signal.h>
#include <inttypes.h>
#include <byteswap.h>
#include <time.h>
#include "s2n-bignum.h"

int usage()
{
    fprintf(stderr, "Usage: wvk offset PUBLIC_KEY PREFIX SKIP LIMIT\n");
    fprintf(stderr, "Find PREFIX offset from PUBLIC_KEY after skipping SKIP steps and making up to LIMIT steps\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Usage: wvk add OFFSET PREFIX\n");
    fprintf(stderr, "Add OFFSET to the private key read from stdin and verify new public key has expected PREFIX\n");
    return 2;
}

static const int batch_size = 4 * 1024; // must be even

typedef uint64_t field_Element[4];

static const int X = 0;
static const int Y = 1;
static const int Z = 2;
static const int T = 3;
static const int XY = 2;

typedef uint64_t PointXY[8];
typedef field_Element Point[4];  // X, Y, Z, T
typedef field_Element affine[3]; // X, Y, XY

void bytes_montgomery(PointXY p, uint8_t *dst);
void vectorDivision(field_Element *x, field_Element *y, field_Element *u, int n);

void PointXY_FromPoint(PointXY v, Point p);
void PointXY_PrintMontgomeryBytes(PointXY p);

void Point_FromXY(Point v, PointXY p);
void Point_PrintMontgomeryBytes(Point p);
void Point_Set(Point v, Point u);

void field_Element_Set(field_Element v, field_Element x);
void field_Element_PrintBytes(field_Element x);

void affine_fromP3(affine v, Point p);
void affine_fromP3zInv(affine v, Point p, field_Element zInv);

uint64_t reverse_bits(uint64_t x);

int base64_encode(const uint8_t *src, int srclen, char *dst);
int base64_decode(const char *src, int srclen, uint8_t *dst);

static field_Element ZERO = {
    UINT64_C(0),
    UINT64_C(0),
    UINT64_C(0),
    UINT64_C(0),
};

static field_Element ONE = {
    UINT64_C(1),
    UINT64_C(0),
    UINT64_C(0),
    UINT64_C(0),
};

static field_Element scalar_offset = {
    UINT64_C(8),
    UINT64_C(0),
    UINT64_C(0),
    UINT64_C(0),
};

static Point point_offset; // b * scalar_offset

static uint8_t prefix_bytes[32];
static int prefix_bits;

void init()
{
    PointXY t;
    edwards25519_scalarmulbase(t, scalar_offset);
    Point_FromXY(point_offset, t);
}

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

    init();
    signal(SIGINT, signal_handler);

    field_Element public_key;
    {
        if (strlen(arg_public_key) == 44)
        {
            char buf[32];
            base64_decode(arg_public_key, 44, buf);
            bignum_fromlebytes_4(public_key, buf);
        }
        else
        {
            fprintf(stderr, "Invalid public key\n");
            return usage();
        }
    }

    uint64_t skipUint = strtoul(arg_skip, NULL, 10);
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

    Point p0;
    {
        PointXY p0xy;
        uint8_t buf32[32];
        field_Element u;
        field_Element_Set(u, public_key);
        // y = (u - 1) / (u + 1)
        field_Element y, n, d, r;

        bignum_sub_p25519(n, u, ONE);
        bignum_add_p25519(d, u, ONE);
        bignum_inv_p25519(r, d);
        bignum_mul_p25519(y, n, r);
        bignum_tolebytes_4(buf32, y);

        if (edwards25519_decode(p0xy, buf32))
        {
            fprintf(stderr, "Invalid public key\n");
            exit(1);
        }
        Point_FromXY(p0, p0xy);
    }

    Point p;
    {
        field_Element skip = {
            skipUint,
            UINT64_C(0),
            UINT64_C(0),
            UINT64_C(0),
        };
        field_Element t;
        bignum_mul_p25519(t, skip, scalar_offset);

        PointXY pt;
        edwards25519_scalarmulbase(pt, t);

        Point skip_offset;
        Point_FromXY(skip_offset, pt);
        edwards25519_epadd((uint64_t *)p, (uint64_t *)p0, (uint64_t *)skip_offset);
    }

    uint64_t n = skipUint;

    field_Element ua[batch_size + 2];
    field_Element ub[batch_size + 2];
    field_Element u[batch_size + 2];

    affine offsets[batch_size / 2];
    Point poi;
    Point_Set(poi, point_offset);
    for (int i = 0; i < batch_size / 2 - 1; i++)
    {
        affine_fromP3(offsets[i], poi);
        edwards25519_epadd((uint64_t *)poi, (uint64_t *)poi, (uint64_t *)point_offset);
    }
    Point batch_offset;
    Point_Set(batch_offset, point_offset);

    affine_fromP3(offsets[batch_size / 2 - 1], poi);
    edwards25519_epadd((uint64_t *)batch_offset, (uint64_t *)batch_offset, (uint64_t *)poi);
    edwards25519_epadd((uint64_t *)batch_offset, (uint64_t *)batch_offset, (uint64_t *)poi);
    edwards25519_epadd((uint64_t *)p, (uint64_t *)p, (uint64_t *)poi);

    affine pa;
    affine_fromP3(pa, p);

    field_Element_Set(ua[batch_size + 1], ONE);

    field_Element num, den, x1y2, y1x2;

    clock_t start_time = clock();
    while (!interrupted)
    {
        for (int i = 0; i < batch_size / 2; i++)
        {
            field_Element *p1 = pa;
            field_Element *p2 = offsets[i];

            bignum_mul_p25519(x1y2, p1[X], p2[Y]);
            bignum_mul_p25519(y1x2, p1[Y], p2[X]);

            bignum_sub_p25519(num, p1[XY], p2[XY]);
            bignum_sub_p25519(den, x1y2, y1x2);
            bignum_add_p25519(ua[batch_size / 2 + 1 + i], den, num);
            bignum_sub_p25519(ub[batch_size / 2 + 1 + i], den, num);

            bignum_add_p25519(num, p1[XY], p2[XY]);
            bignum_add_p25519(den, x1y2, y1x2);
            bignum_add_p25519(ua[batch_size / 2 - 1 - i], den, num);
            bignum_sub_p25519(ub[batch_size / 2 - 1 - i], den, num);
        }
        bignum_add_p25519(ua[batch_size / 2], ONE, pa[Y]);
        bignum_sub_p25519(ub[batch_size / 2], ONE, pa[Y]);

        edwards25519_epadd((uint64_t *)p, (uint64_t *)p, (uint64_t *)batch_offset);

        field_Element_Set(ub[batch_size + 1], p[Z]);
        uint64_t *pZ_inv = u[batch_size + 1];

        vectorDivision(ua, ub, u, batch_size + 2);

        for (int i = 0; i < batch_size + 1; i++)
        {
            if ((u[i][0] & mask) == prefix_match)
            {
                n += i;
                goto found;
            }
        }

        n += batch_size + 1;
        affine_fromP3zInv(pa, p, pZ_inv);

        if (limit > 0)
        {
            if (limit <= batch_size + 1)
            {
                break;
            }
            limit -= batch_size + 1;
        }
    }

found:
    double seconds = (double)(clock() - start_time) / (double)CLOCKS_PER_SEC;

    printf("%lu\n", n);

    fprintf(stderr, "seconds: %0.0f\n", seconds);
    fprintf(stderr, "attempts/s: %0.0f\n", (n - skipUint) / seconds);

    if (interrupted)
    {
        exit(3);
    }
}

// add OFFSET PREFIX
int cmd_add(int argc, char *argv[])
{
    const char *arg_offset = argv[1];
    const char *arg_prefix = argv[2];

    field_Element s0;
    {
        uint8_t buf32[32];
        char buf[45];
        if (fgets(buf, 45, stdin) == NULL)
        {
            fprintf(stderr, "Failed to read private key from stdin\n");
            usage();
        }
        buf[44] = '\0';
        base64_decode(buf, 44, buf32);
        bignum_fromlebytes_4(s0, buf32);
    }

    uint64_t offset = strtoul(arg_offset, NULL, 10);

    field_Element so;
    field_Element oo = {
        offset,
        UINT64_C(0),
        UINT64_C(0),
        UINT64_C(0),
    };
    bignum_mul_p25519(so, scalar_offset, oo);

    field_Element sp, sm;
    bignum_add_p25519(sp, s0, so);
    bignum_sub_p25519(sm, s0, so);

    PointXY pp, pm;
    edwards25519_scalarmulbase(pp, sp);
    edwards25519_scalarmulbase(pm, sm);

    uint8_t buf32[32];
    uint8_t buf45[45];
    buf45[44] = 0;

    bytes_montgomery(pp, buf32);
    base64_encode(buf32, 32, buf45);
    if (!strncmp(buf45, arg_prefix, strlen(arg_prefix)))
    {
        // key bytes
        field_Element_PrintBytes(sp);
        return 0;
    }

    bytes_montgomery(pm, buf32);
    base64_encode(buf32, 32, buf45);
    if (!strncmp(buf45, arg_prefix, strlen(arg_prefix)))
    {
        field_Element_PrintBytes(sm);
        return 0;
    }

    fprintf(stderr, "Prefix mismatch\n");
    return 1;
}

void bytes_montgomery(PointXY p, uint8_t *dst)
{
    // u = (1 + y) / (1 - y)
    uint64_t *y = p + 4;
    field_Element n, d, t, u;

    bignum_add_p25519(n, ONE, y);
    bignum_sub_p25519(d, ONE, y);
    bignum_inv_p25519(t, d);
    bignum_mul_p25519(u, n, t);

    bignum_tolebytes_4(dst, u);
}

// vectorDivision calculates u = x / y
//
// It uses:
//
//	4*(n-1)+1 multiplications
//	1 invert = ~265 multiplications
//
// Complexity: 262M + 4M*n
//
// Simultaneous field divisions: an extension of Montgomery's trick
// David G. Harris
// https://eprint.iacr.org/2008/199.pdf
void vectorDivision(field_Element *x, field_Element *y, field_Element *u, int n)
{
    field_Element a, b;
    uint64_t *ri = a, *ri_1 = b; // r[i], r[i-1]

    field_Element_Set(ri_1, y[0]);
    for (int i = 1; i < n; i++)
    {
        // TODO: For reasons to be explored
        // the following two multiplications in that exact order and
        // use of swapping pointers here and in the second loop
        // achieve better performance than a straightforward implementation:
        //
        //     field_Element py; // y[0]*y[1]*...*y[n]
        //     field_Element_Set(py, y[0]);
        //     for (int i = 1; i < n; i++)
        //     {
        //         bignum_mul_p25519(u[i], py, x[i]);
        //         bignum_mul_p25519(py, py, y[i]);
        //     }
        //
        bignum_mul_p25519(ri, ri_1, y[i]);
        bignum_mul_p25519(u[i], ri_1, x[i]);

        uint64_t *t = ri;
        ri = ri_1;
        ri_1 = t;
    }

    field_Element I, J;
    uint64_t *pI = I, *pJ = J;
    bignum_inv_p25519(pI, ri_1); // r[n-1]

    for (int i = n - 1; i > 0; i--)
    {
        bignum_mul_p25519(pJ, pI, y[i]);
        bignum_mul_p25519(u[i], pI, u[i]);

        uint64_t *t = pJ;
        pJ = pI;
        pI = t;
    }
    bignum_mul_p25519(u[0], pI, x[0]);
}

void Point_PrintMontgomeryBytes(Point p)
{
    field_Element n, d, t, u;
    // u = (1 + y) / (1 - y) = (Z + Y) / (Z - Y)
    bignum_add_p25519(n, p[Z], p[Y]);
    bignum_sub_p25519(d, p[Z], p[Y]);
    bignum_inv_p25519(t, d);
    bignum_mul_p25519(u, n, t);

    field_Element_PrintBytes(u);
}

void Point_Set(Point v, Point u)
{
    field_Element_Set(v[X], u[X]);
    field_Element_Set(v[Y], u[Y]);
    field_Element_Set(v[Z], u[Z]);
    field_Element_Set(v[T], u[T]);
}

void Point_FromXY(Point v, PointXY p)
{
    field_Element_Set(v[X], p);
    field_Element_Set(v[Y], p + 4);
    field_Element_Set(v[Z], ONE);
    bignum_mul_p25519(v[T], v[X], v[Y]);
}

void PointXY_FromPoint(PointXY v, Point p)
{
    field_Element t;
    bignum_inv_p25519(t, p[Z]);        // t = 1 / Z
    bignum_mul_p25519(v, p[X], t);     // x = X / Z
    bignum_mul_p25519(v + 4, p[Y], t); // y = Y / Z
}

void PointXY_PrintMontgomeryBytes(PointXY p)
{
    uint8_t buf32[32];
    uint8_t buf45[45];
    buf45[44] = 0;

    bytes_montgomery(p, buf32);
    base64_encode(buf32, 32, buf45);
    printf("%s\n", buf45);
}

void field_Element_Set(field_Element v, field_Element x)
{
    v[0] = x[0];
    v[1] = x[1];
    v[2] = x[2];
    v[3] = x[3];
}

void field_Element_PrintBytes(field_Element x)
{
    uint8_t buf32[32];
    uint8_t buf45[45];
    buf45[44] = 0;

    bignum_tolebytes_4(buf32, x);
    base64_encode(buf32, 32, buf45);
    printf("%s\n", buf45);
}

void affine_fromP3(affine v, Point p)
{
    field_Element zInv;
    bignum_inv_p25519(zInv, p[Z]);
    bignum_mul_p25519(v[X], p[X], zInv);
    bignum_mul_p25519(v[Y], p[Y], zInv);
    bignum_mul_p25519(v[XY], v[X], v[Y]);
}

void affine_fromP3zInv(affine v, Point p, field_Element zInv)
{
    bignum_mul_p25519(v[X], p[X], zInv);
    bignum_mul_p25519(v[Y], p[Y], zInv);
    bignum_mul_p25519(v[XY], v[X], v[Y]);
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

static const char *base64_digits = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/\0";

int base64_decode(const char *src, int srclen, uint8_t *dst)
{
    int j = 0;
    int shift = 0;
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
        dst[j++] = base64_digits[((src[i] << 4) & 0b110000) | ((src[++i] >> 4) & 0b001111)];
        dst[j++] = base64_digits[(src[i] << 2) & 0b111111];
        dst[j++] = '=';
        break;
    default:
        return -1;
    }
    return j;
}

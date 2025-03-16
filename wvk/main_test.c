#include "main.c"

#define expect(expr) ((expr) ? (void)(0) : fail(#expr, __FILE__, __LINE__, __func__))

void fail(const char *assertion, const char *file, unsigned int line, const char *function)
{
    fprintf(stderr, "%s:%d: %s failed in %s\n", file, line, assertion, function);
    exit(1);
}

void test_base64_decode()
{
    uint8_t buf[32];

    expect(base64_decode("B", 1, buf) == 1);
    expect(buf[0] == 0b00000100);

    expect(base64_decode("BB", 2, buf) == 2);
    expect(buf[0] == 0b00000100);
    expect(buf[1] == 0b00010000);

    expect(base64_decode("BBB", 3, buf) == 3);
    expect(buf[0] == 0b00000100);
    expect(buf[1] == 0b00010000);
    expect(buf[2] == 0b01000000);
}

void test_base64_encode()
{
    uint8_t buf[45];

    expect(base64_encode("M", 1, buf) == 4);
    expect(strncmp(buf, "TQ==", 4) == 0);

    expect(base64_encode("Ma", 2, buf) == 4);
    expect(strncmp(buf, "TWE=", 4) == 0);

    expect(base64_encode("\x01\x01\x01", 3, buf) == 4);
    expect(strncmp(buf, "AQEB", 4) == 0);

    expect(base64_encode("\x01\x01\x01\x01", 4, buf) == 8);
    expect(strncmp(buf, "AQEBAQ==", 8) == 0);

    const char *key = "YJbwTLGyZxJCuVETzG7VQb1DGqsKLGAsa07mEoERTFM=";
    expect(base64_decode(key, 44, buf) == 33);

    uint8_t buf45[45];
    base64_encode(buf, 32, buf45);
    buf45[44] = 0;

    expect(!strncmp(buf45, key, 44));
}

int test_all(int argc, char *argv[])
{
    test_base64_decode();
    test_base64_encode();

    exit(0);
}

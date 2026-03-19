/*
 * ML-DSA-65 test program.
 * - NIST ACVP keygen test vector (tcId 26, ML-DSA-65)
 * - Deterministic sign + verify round-trip
 */

#include "ml-dsa-65.h"
#include "sha3.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

/* ---- Hex utilities ---- */

static int hex2bin(uint8_t *out, const char *hex, size_t out_len) {
  for (size_t i = 0; i < out_len; i++) {
    unsigned int hi, lo;
    if (sscanf(hex + 2 * i, "%1x", &hi) != 1) return -1;
    if (sscanf(hex + 2 * i + 1, "%1x", &lo) != 1) return -1;
    out[i] = (uint8_t)((hi << 4) | lo);
  }
  return 0;
}

static void print_hex(const char *label, const uint8_t *data, size_t len) {
  printf("%s: ", label);
  for (size_t i = 0; i < len; i++)
    printf("%02x", data[i]);
  printf("\n");
}

/* ---- Test 1: KeyGen + Sign + Verify round-trip ---- */

static int test_roundtrip(void) {
  /* ACVP ML-DSA-65 keyGen test vector, tcId 26 */
  const char *seed_hex =
    "1BD67DC782B2958E189E315C040DD1F64C8AB232A6A170E1A7A52C33F10851B1";

  uint8_t seed[32];
  hex2bin(seed, seed_hex, 32);

  static uint8_t pk[MLDSA_PK_BYTES];
  static uint8_t sk[MLDSA_SK_BYTES];
  static uint8_t sig[MLDSA_SIG_BYTES];
  size_t sig_len = 0;

  printf("=== Test 1: KeyGen + Sign + Verify round-trip ===\n");

  /* KeyGen */
  printf("Running KeyGen...\n");
  int rc = ml_dsa_65_keygen(pk, sk, seed);
  if (rc != 0) {
    printf("FAIL: KeyGen returned %d\n", rc);
    return 1;
  }

  /* Hash pk and sk for KAT comparison */
  {
    uint8_t pk_hash[32], sk_hash[32];
    sha3_256_raw(pk, MLDSA_PK_BYTES, pk_hash);
    sha3_256_raw(sk, MLDSA_SK_BYTES, sk_hash);
    print_hex("  pk SHA3-256", pk_hash, 32);
    print_hex("  sk SHA3-256", sk_hash, 32);
  }

  /* Sign */
  const uint8_t msg[] = "test message for ML-DSA-65";
  size_t msg_len = sizeof(msg) - 1;

  printf("Running Sign (deterministic)...\n");
  rc = ml_dsa_65_sign(sig, &sig_len, msg, msg_len, NULL, 0, sk);
  if (rc != 0) {
    printf("FAIL: Sign returned %d\n", rc);
    return 1;
  }
  printf("  Signature length: %zu bytes (expected %d)\n",
         sig_len, MLDSA_SIG_BYTES);

  {
    uint8_t sig_hash[32];
    sha3_256_raw(sig, sig_len, sig_hash);
    print_hex("  sig SHA3-256", sig_hash, 32);
  }

  /* Verify */
  printf("Running Verify...\n");
  rc = ml_dsa_65_verify(msg, msg_len, sig, sig_len, NULL, 0, pk);
  if (rc != 0) {
    printf("FAIL: Verify returned %d\n", rc);
    return 1;
  }
  printf("  Verify: PASS\n");

  /* Verify with wrong message should fail */
  const uint8_t bad_msg[] = "wrong message";
  rc = ml_dsa_65_verify(bad_msg, sizeof(bad_msg) - 1, sig, sig_len, NULL, 0, pk);
  if (rc == 0) {
    printf("FAIL: Verify should reject wrong message\n");
    return 1;
  }
  printf("  Reject bad msg: PASS\n");

  /* Verify with corrupted signature should fail */
  sig[100] ^= 0x01;
  rc = ml_dsa_65_verify(msg, msg_len, sig, sig_len, NULL, 0, pk);
  if (rc == 0) {
    printf("FAIL: Verify should reject corrupted signature\n");
    return 1;
  }
  printf("  Reject bad sig: PASS\n");
  sig[100] ^= 0x01; /* restore */

  return 0;
}

/* ---- Test 2: Deterministic signing is reproducible ---- */

static int test_deterministic(void) {
  const char *seed_hex =
    "CA464DD4C09BA7346057527285D84AAC437DB1525EF72403D93D8E0E9301E9A0";

  uint8_t seed[32];
  hex2bin(seed, seed_hex, 32);

  static uint8_t pk[MLDSA_PK_BYTES];
  static uint8_t sk[MLDSA_SK_BYTES];
  static uint8_t sig1[MLDSA_SIG_BYTES];
  static uint8_t sig2[MLDSA_SIG_BYTES];
  size_t sig_len1 = 0, sig_len2 = 0;

  printf("\n=== Test 2: Deterministic signing reproducibility ===\n");

  ml_dsa_65_keygen(pk, sk, seed);

  const uint8_t msg[] = "determinism check";
  size_t msg_len = sizeof(msg) - 1;

  ml_dsa_65_sign(sig1, &sig_len1, msg, msg_len, NULL, 0, sk);
  ml_dsa_65_sign(sig2, &sig_len2, msg, msg_len, NULL, 0, sk);

  if (sig_len1 != sig_len2 || memcmp(sig1, sig2, sig_len1) != 0) {
    printf("FAIL: Deterministic signatures differ\n");
    return 1;
  }
  printf("  Deterministic: PASS\n");

  /* Verify both */
  int rc = ml_dsa_65_verify(msg, msg_len, sig1, sig_len1, NULL, 0, pk);
  if (rc != 0) {
    printf("FAIL: Verify returned %d\n", rc);
    return 1;
  }
  printf("  Verify: PASS\n");
  return 0;
}

/* ---- Test 3: Context string support ---- */

static int test_context_string(void) {
  const char *seed_hex =
    "684201C617E77778DBC6F0634E336C275C27401248440E2B0D01846746FB1CD0";

  uint8_t seed[32];
  hex2bin(seed, seed_hex, 32);

  static uint8_t pk[MLDSA_PK_BYTES];
  static uint8_t sk[MLDSA_SK_BYTES];
  static uint8_t sig[MLDSA_SIG_BYTES];
  size_t sig_len = 0;

  printf("\n=== Test 3: Context string support ===\n");

  ml_dsa_65_keygen(pk, sk, seed);

  const uint8_t msg[] = "context test message";
  size_t msg_len = sizeof(msg) - 1;
  const uint8_t ctx[] = "my-app-v1";
  size_t ctx_len = sizeof(ctx) - 1;

  int rc = ml_dsa_65_sign(sig, &sig_len, msg, msg_len, ctx, ctx_len, sk);
  if (rc != 0) {
    printf("FAIL: Sign with context returned %d\n", rc);
    return 1;
  }

  rc = ml_dsa_65_verify(msg, msg_len, sig, sig_len, ctx, ctx_len, pk);
  if (rc != 0) {
    printf("FAIL: Verify with context returned %d\n", rc);
    return 1;
  }
  printf("  Sign+Verify with context: PASS\n");

  /* Verify with wrong context should fail */
  const uint8_t bad_ctx[] = "wrong-ctx";
  rc = ml_dsa_65_verify(msg, msg_len, sig, sig_len, bad_ctx, sizeof(bad_ctx) - 1, pk);
  if (rc == 0) {
    printf("FAIL: Verify should reject wrong context\n");
    return 1;
  }
  printf("  Reject wrong context: PASS\n");

  /* Verify with empty context should fail */
  rc = ml_dsa_65_verify(msg, msg_len, sig, sig_len, NULL, 0, pk);
  if (rc == 0) {
    printf("FAIL: Verify should reject empty context\n");
    return 1;
  }
  printf("  Reject empty context: PASS\n");

  return 0;
}

/* ---- Main ---- */

int main(void) {
  int failures = 0;

  failures += test_roundtrip();
  failures += test_deterministic();
  failures += test_context_string();

  printf("\n=== Summary: %s (%d failure%s) ===\n",
         failures == 0 ? "ALL PASSED" : "SOME FAILED",
         failures, failures == 1 ? "" : "s");

  return failures;
}

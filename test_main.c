/*
 * ML-DSA-65 test program.
 * - NIST ACVP keygen test vector (tcId 26, ML-DSA-65)
 * - Deterministic sign + verify round-trip
 */

#include "ml-dsa-65.h"
#define MLD_CONFIG_API_PARAMETER_SET 65
#define MLD_CONFIG_API_NAMESPACE_PREFIX mldsa
#include "mldsa_native.h"
#include "sha3.h"
#include "test_vectors.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

typedef struct {
  const char *seed_hex;
  const char *msg;
  const char *ctx;
} ref_case_t;

static const ref_case_t ref_cases[] = {
  { ACVP_TC26_SEED_HEX, "cross-check signature 0", "interop-ctx-0" },
  { "B850D898A3D3D11C4E64ADE5A86FFED951B237C60D2A67A2DEF0A792B8F6990D", "cross-check signature 1", "interop-ctx-1" },
  { "455ECBD3C4A9EFB75A302DF08E770BF79E8605DC13ED57D7319AA6BFD1B6496B", "cross-check signature 2", "interop-ctx-2" },
  { "DDC3DE6AAA57CCF19272FB4CC76D933D292D11921CA93F4AB3DBE18AFD9A5DF0", "cross-check signature 3", "interop-ctx-3" },
  { "CA464DD4C09BA7346057527285D84AAC437DB1525EF72403D93D8E0E9301E9A0", "cross-check signature 4", "interop-ctx-4" },
  { "684201C617E77778DBC6F0634E336C275C27401248440E2B0D01846746FB1CD0", "cross-check signature 5", "interop-ctx-5" },
  { "2DFC84C9589EA2C45124288B86CDE151797FA6F0C94EB762501381E9CEF707C3", "cross-check signature 6", "interop-ctx-6" },
  { "63F66260BCF8CD12F12745C3A45E29BDCFB3B5789CBDD7EA37D96B3155DB4997", "cross-check signature 7", "interop-ctx-7" }
};

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

static int expect_buf_eq(const char *label,
                         const uint8_t *got,
                         const uint8_t *expected,
                         size_t len) {
  if (memcmp(got, expected, len) == 0)
    return 0;

  printf("FAIL: %s mismatch\n", label);
  return 1;
}

static int ref_keygen(uint8_t *pk, uint8_t *sk, const uint8_t *seed) {
  return MLD_API_NAMESPACE(keypair_internal)(pk, sk, seed);
}

static int ref_sign_deterministic(uint8_t *sig, size_t *sig_len,
                                  const uint8_t *msg, size_t msg_len,
                                  const uint8_t *ctx, size_t ctx_len,
                                  const uint8_t *sk) {
  uint8_t prefix[257];
  uint8_t rnd[32] = {0};
  size_t prefix_len;

  if (ctx_len > 255)
    return -1;

  prefix[0] = 0x00;
  prefix[1] = (uint8_t)ctx_len;
  if (ctx_len > 0 && ctx != NULL)
    memcpy(prefix + 2, ctx, ctx_len);
  prefix_len = 2 + ctx_len;

  if (prefix_len == 0)
    return -1;

  return MLD_API_NAMESPACE(signature_internal)(
      sig, sig_len, msg, msg_len, prefix, prefix_len, rnd, sk, 0);
}

static int ref_verify_sig(const uint8_t *sig, size_t sig_len,
                          const uint8_t *msg, size_t msg_len,
                          const uint8_t *ctx, size_t ctx_len,
                          const uint8_t *pk) {
  return MLD_API_NAMESPACE(verify)(sig, sig_len, msg, msg_len, ctx, ctx_len, pk);
}

static int test_acvp_keygen_vector(void) {
  uint8_t seed[32];
  static uint8_t pk[MLDSA_PK_BYTES];
  static uint8_t sk[MLDSA_SK_BYTES];
  static uint8_t expected_pk[MLDSA_PK_BYTES];
  static uint8_t expected_sk[MLDSA_SK_BYTES];

  printf("=== Test 1: ACVP keyGen tcId 26 exact vector ===\n");

  if (hex2bin(seed, ACVP_TC26_SEED_HEX, sizeof(seed)) != 0)
    return 1;
  if (hex2bin(expected_pk, ACVP_TC26_PK_HEX, sizeof(expected_pk)) != 0)
    return 1;
  if (hex2bin(expected_sk, ACVP_TC26_SK_HEX, sizeof(expected_sk)) != 0)
    return 1;

  if (ml_dsa_65_keygen(pk, sk, NULL, seed) != 0) {
    printf("FAIL: KeyGen returned non-zero\n");
    return 1;
  }

  if (expect_buf_eq("ACVP pk", pk, expected_pk, sizeof(pk)) != 0)
    return 1;
  if (expect_buf_eq("ACVP sk", sk, expected_sk, sizeof(sk)) != 0)
    return 1;

  printf("  ACVP keyGen exact match: PASS\n");
  return 0;
}

static int test_reference_keygen_equivalence(void) {
  printf("\n=== Test 2: Reference keygen equivalence ===\n");

  for (size_t i = 0; i < sizeof(ref_cases) / sizeof(ref_cases[0]); i++) {
    uint8_t seed[32];
    static uint8_t pk[MLDSA_PK_BYTES];
    static uint8_t sk[MLDSA_SK_BYTES];
    static uint8_t ref_pk[MLDSA_PK_BYTES];
    static uint8_t ref_sk[MLDSA_SK_BYTES];

    if (hex2bin(seed, ref_cases[i].seed_hex, sizeof(seed)) != 0)
      return 1;
    if (ml_dsa_65_keygen(pk, sk, NULL, seed) != 0)
      return 1;
    if (ref_keygen(ref_pk, ref_sk, seed) != 0) {
      printf("FAIL: Reference keyGen returned non-zero for case %zu\n", i + 1);
      return 1;
    }
    if (expect_buf_eq("reference pk", pk, ref_pk, sizeof(pk)) != 0)
      return 1;
    if (expect_buf_eq("reference sk", sk, ref_sk, sizeof(sk)) != 0)
      return 1;
  }

  printf("  Reference keygen: PASS\n");
  return 0;
}

/* ---- Test 1: KeyGen + Sign + Verify round-trip ---- */

static int test_roundtrip(void) {
  const char *seed_hex =
    ACVP_TC26_SEED_HEX;

  uint8_t seed[32];
  hex2bin(seed, seed_hex, 32);

  static uint8_t pk[MLDSA_PK_BYTES];
  static uint8_t sk[MLDSA_SK_BYTES];
  static uint8_t sig[MLDSA_SIG_BYTES];
  size_t sig_len = 0;

  printf("=== Test 1: KeyGen + Sign + Verify round-trip ===\n");

  /* KeyGen */
  printf("Running KeyGen...\n");
  int rc = ml_dsa_65_keygen(pk, sk, NULL, seed);
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

static int test_reference_sign_interop(void) {
  printf("\n=== Test 4: Reference sign/verify interop ===\n");

  for (size_t i = 0; i < sizeof(ref_cases) / sizeof(ref_cases[0]); i++) {
    uint8_t seed[32];
    static uint8_t pk[MLDSA_PK_BYTES];
    static uint8_t sk[MLDSA_SK_BYTES];
    static uint8_t ref_pk[MLDSA_PK_BYTES];
    static uint8_t ref_sk[MLDSA_SK_BYTES];
    static uint8_t sig[MLDSA_SIG_BYTES];
    static uint8_t ref_sig[MLDSA_SIG_BYTES];
    size_t sig_len = 0;
    size_t ref_sig_len = 0;
    const uint8_t *msg = (const uint8_t *)ref_cases[i].msg;
    size_t msg_len = strlen(ref_cases[i].msg);
    const uint8_t *ctx = (const uint8_t *)ref_cases[i].ctx;
    size_t ctx_len = strlen(ref_cases[i].ctx);

    if (hex2bin(seed, ref_cases[i].seed_hex, 32) != 0)
      return 1;
    if (ml_dsa_65_keygen(pk, sk, NULL, seed) != 0)
      return 1;
    if (ref_keygen(ref_pk, ref_sk, seed) != 0)
      return 1;

    if (ml_dsa_65_sign(sig, &sig_len, msg, msg_len, ctx, ctx_len, sk) != 0) {
      printf("FAIL: local sign returned non-zero\n");
      return 1;
    }
    if (ref_sign_deterministic(ref_sig, &ref_sig_len, msg, msg_len, ctx, ctx_len, ref_sk) != 0) {
      printf("FAIL: reference sign returned non-zero\n");
      return 1;
    }
    if (sig_len != ref_sig_len) {
      printf("FAIL: signature length mismatch\n");
      return 1;
    }
    if (expect_buf_eq("reference signature", sig, ref_sig, sig_len) != 0)
      return 1;
    if (ml_dsa_65_verify(msg, msg_len, ref_sig, ref_sig_len, ctx, ctx_len, pk) != 0) {
      printf("FAIL: local verify rejected reference signature\n");
      return 1;
    }
    if (ref_verify_sig(sig, sig_len, msg, msg_len, ctx, ctx_len, ref_pk) != 0) {
      printf("FAIL: reference verify rejected local signature\n");
      return 1;
    }

    {
      const uint8_t bad_msg[] = "cross-check signature mismatch";
      if (ml_dsa_65_verify(bad_msg, sizeof(bad_msg) - 1, ref_sig, ref_sig_len, ctx, ctx_len, pk) == 0) {
        printf("FAIL: local verify accepted mismatched message\n");
        return 1;
      }
      if (ref_verify_sig(sig, sig_len, bad_msg, sizeof(bad_msg) - 1, ctx, ctx_len, ref_pk) == 0) {
        printf("FAIL: reference verify accepted mismatched message\n");
        return 1;
      }
    }
  }

  printf("  Reference sign/verify interop: PASS\n");
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

  ml_dsa_65_keygen(pk, sk, NULL, seed);

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

  ml_dsa_65_keygen(pk, sk, NULL, seed);

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

/* ---- Test 5: Randomized sign/verify interop ---- */

static int test_random_interop(void) {
  int nrounds = 32;
  static uint8_t pk[MLDSA_PK_BYTES];
  static uint8_t sk[MLDSA_SK_BYTES];
  static uint8_t ref_pk[MLDSA_PK_BYTES];
  static uint8_t ref_sk[MLDSA_SK_BYTES];
  static uint8_t sig[MLDSA_SIG_BYTES];
  static uint8_t ref_sig[MLDSA_SIG_BYTES];
  uint8_t seed[32];
  uint8_t msg[128];
  uint8_t ctx[16];

  printf("\n=== Test 5: Randomized sign/verify interop (%d rounds) ===\n",
         nrounds);

  /* Use SHAKE256 as a deterministic PRNG seeded from a fixed value,
   * so the test is reproducible but covers many inputs. */
  SHA3_CTX_T prng;
  uint8_t prng_seed[32] = "pqc-test-random-interop-seed!!!!";
  shake256_init(&prng);
  shake_update(&prng, prng_seed, 32);
  shake_finalize(&prng);

  for (int round = 0; round < nrounds; round++) {
    /* Generate random seed, message, context from PRNG */
    shake_squeeze(&prng, seed, 32);
    uint8_t lens[2];
    shake_squeeze(&prng, lens, 2);
    size_t msg_len = 1 + (lens[0] % sizeof(msg));  /* 1..128 */
    size_t ctx_len = lens[1] % (sizeof(ctx) + 1);  /* 0..16 */
    shake_squeeze(&prng, msg, msg_len);
    if (ctx_len > 0)
      shake_squeeze(&prng, ctx, ctx_len);

    /* KeyGen: both implementations */
    if (ml_dsa_65_keygen(pk, sk, NULL, seed) != 0) {
      printf("FAIL round %d: local keygen\n", round);
      return 1;
    }
    if (ref_keygen(ref_pk, ref_sk, seed) != 0) {
      printf("FAIL round %d: ref keygen\n", round);
      return 1;
    }
    if (memcmp(pk, ref_pk, MLDSA_PK_BYTES) != 0 ||
        memcmp(sk, ref_sk, MLDSA_SK_BYTES) != 0) {
      printf("FAIL round %d: keygen mismatch\n", round);
      return 1;
    }

    /* Sign: both implementations, compare signatures */
    size_t sig_len = 0, ref_sig_len = 0;
    if (ml_dsa_65_sign(sig, &sig_len, msg, msg_len,
                       ctx_len > 0 ? ctx : NULL, ctx_len, sk) != 0) {
      printf("FAIL round %d: local sign\n", round);
      return 1;
    }
    if (ref_sign_deterministic(ref_sig, &ref_sig_len, msg, msg_len,
                               ctx_len > 0 ? ctx : NULL, ctx_len,
                               ref_sk) != 0) {
      printf("FAIL round %d: ref sign\n", round);
      return 1;
    }
    if (sig_len != ref_sig_len || memcmp(sig, ref_sig, sig_len) != 0) {
      printf("FAIL round %d: signature mismatch\n", round);
      return 1;
    }

    /* Cross-verify */
    if (ml_dsa_65_verify(msg, msg_len, ref_sig, ref_sig_len,
                         ctx_len > 0 ? ctx : NULL, ctx_len, pk) != 0) {
      printf("FAIL round %d: local verify ref sig\n", round);
      return 1;
    }
    if (ref_verify_sig(sig, sig_len, msg, msg_len,
                       ctx_len > 0 ? ctx : NULL, ctx_len, ref_pk) != 0) {
      printf("FAIL round %d: ref verify local sig\n", round);
      return 1;
    }
  }

  printf("  Randomized interop (%d rounds): PASS\n", nrounds);
  return 0;
}

/* ---- Main ---- */

int main(void) {
  int failures = 0;

  failures += test_acvp_keygen_vector();
  failures += test_reference_keygen_equivalence();
  failures += test_roundtrip();
  failures += test_deterministic();
  failures += test_reference_sign_interop();
  failures += test_context_string();
  failures += test_random_interop();

  printf("\n=== Summary: %s (%d failure%s) ===\n",
         failures == 0 ? "ALL PASSED" : "SOME FAILED",
         failures, failures == 1 ? "" : "s");

  return failures;
}

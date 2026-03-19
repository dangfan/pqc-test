/*
 * Debug: test NTT round-trip and sign with progress output.
 * We instrument by directly calling keygen then sign and printing intermediate state.
 */
#include "ml-dsa-65.h"
#include "sha3.h"
#include <stdio.h>
#include <string.h>
#include <time.h>

int main(void) {
  static uint8_t pk[MLDSA_PK_BYTES];
  static uint8_t sk[MLDSA_SK_BYTES];
  static uint8_t sig[MLDSA_SIG_BYTES];
  size_t sig_len = 0;

  uint8_t seed[32] = {0};
  seed[0] = 0x42;

  printf("KeyGen...\n");
  clock_t t0 = clock();
  int rc = ml_dsa_65_keygen(pk, sk, seed);
  clock_t t1 = clock();
  printf("KeyGen done: rc=%d, %.3f sec\n", rc,
         (double)(t1 - t0) / CLOCKS_PER_SEC);

  const uint8_t msg[] = "hello";
  printf("Sign...\n");
  t0 = clock();
  rc = ml_dsa_65_sign(sig, &sig_len, msg, 5, NULL, 0, sk);
  t1 = clock();
  printf("Sign done: rc=%d, sig_len=%zu, %.3f sec\n", rc, sig_len,
         (double)(t1 - t0) / CLOCKS_PER_SEC);

  if (rc == 0) {
    printf("Verify...\n");
    rc = ml_dsa_65_verify(msg, 5, sig, sig_len, NULL, 0, pk);
    printf("Verify: rc=%d\n", rc);
  }

  return 0;
}

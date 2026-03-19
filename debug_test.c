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
  printf("Running ml_dsa_65_selftest()...\n");
  clock_t t0 = clock();
  int rc = ml_dsa_65_selftest();
  clock_t t1 = clock();
  printf("selftest: rc=%d (%.3f sec)\n", rc,
         (double)(t1 - t0) / CLOCKS_PER_SEC);
  if (rc == 0)
    printf("PASS\n");
  else
    printf("FAIL (error code %d)\n", rc);
  return rc != 0;
}

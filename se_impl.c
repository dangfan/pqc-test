#include "se.h"
#include <assert.h>
#include <mbedtls/bignum.h>
#include <string.h>

// count of uint32_t data
static int pkeLen;
static uint8_t reg[48][64];
static mbedtls_mpi Rinv;

int rsaPKESetLen(int len) {
  assert(len % 2 == 0);
  pkeLen = len;
  assert(pkeLen <= 64);
  return 0;
}

int eccCalMc64(uint32_t *mc, const uint32_t *p) { return 0; }

int rsaPKEWriteMc(const uint32_t *buf) {
  mbedtls_mpi a, n;
  mbedtls_mpi_init(&a);
  mbedtls_mpi_init(&n);
  mbedtls_mpi_init(&Rinv);

  mbedtls_mpi_read_binary(&n, reg[0], pkeLen * 4);

  uint8_t R[pkeLen * 4 + 1]; // R
  memset(R, 0, pkeLen * 4 + 1);
  R[0] = 1;
  mbedtls_mpi_read_binary(&a, R, pkeLen * 4 + 1); // a = R

  mbedtls_mpi_inv_mod(&Rinv, &a, &n);   // Rinv = a^-1 = R^-1

  mbedtls_mpi_free(&a);
  mbedtls_mpi_free(&n);
  return 0;
}

int eccPKEWriteBuf(int idx, const uint32_t *buf, int len) {
  memcpy(reg[idx], buf, len * 4);
  return 0;
}

int eccPKEReadBuf(uint32_t *buf, int idx, int len) {
  memcpy(buf, reg[idx], len * 4);
  return 0;
}

int sm2MonAdd(int r, int idxA, int idxB) {
  mbedtls_mpi re, a, b, n;
  mbedtls_mpi_init(&re);
  mbedtls_mpi_init(&a);
  mbedtls_mpi_init(&b);
  mbedtls_mpi_init(&n);

  mbedtls_mpi_read_binary(&n, reg[0], pkeLen * 4);
  mbedtls_mpi_read_binary(&a, reg[idxA], pkeLen * 4);
  mbedtls_mpi_read_binary(&b, reg[idxB], pkeLen * 4);

  mbedtls_mpi_add_mpi(&re, &a, &b);
  mbedtls_mpi_mod_mpi(&re, &re, &n);

  mbedtls_mpi_write_binary(&re, reg[r], pkeLen * 4);

  mbedtls_mpi_free(&re);
  mbedtls_mpi_free(&a);
  mbedtls_mpi_free(&b);
  mbedtls_mpi_free(&n);

  return 0;
}

int sm2MonSub(int r, int idxA, int idxB) {
  mbedtls_mpi re, a, b, n;
  mbedtls_mpi_init(&re);
  mbedtls_mpi_init(&a);
  mbedtls_mpi_init(&b);
  mbedtls_mpi_init(&n);

  mbedtls_mpi_read_binary(&n, reg[0], pkeLen * 4);
  mbedtls_mpi_read_binary(&a, reg[idxA], pkeLen * 4);
  mbedtls_mpi_read_binary(&b, reg[idxB], pkeLen * 4);

  mbedtls_mpi_sub_mpi(&re, &a, &b);
  mbedtls_mpi_mod_mpi(&re, &re, &n);

  mbedtls_mpi_write_binary(&re, reg[r], pkeLen * 4);

  mbedtls_mpi_free(&re);
  mbedtls_mpi_free(&a);
  mbedtls_mpi_free(&b);
  mbedtls_mpi_free(&n);

  return 0;
}

int sm2MonMul(int r, int idxA, int idxB) {
  mbedtls_mpi re, a, b, n;
  mbedtls_mpi_init(&re);
  mbedtls_mpi_init(&a);
  mbedtls_mpi_init(&b);
  mbedtls_mpi_init(&n);

  mbedtls_mpi_read_binary(&n, reg[0], pkeLen * 4);
  mbedtls_mpi_read_binary(&a, reg[idxA], pkeLen * 4);
  mbedtls_mpi_read_binary(&b, reg[idxB], pkeLen * 4);

  mbedtls_mpi_mul_mpi(&re, &a, &b); // re = a * b
  mbedtls_mpi_mul_mpi(&re, &re, &Rinv); // re = re * Rinv
  mbedtls_mpi_mod_mpi(&re, &re, &n);

  mbedtls_mpi_write_binary(&re, reg[r], pkeLen * 4);

  mbedtls_mpi_free(&re);
  mbedtls_mpi_free(&a);
  mbedtls_mpi_free(&b);
  mbedtls_mpi_free(&n);

  return 0;
}

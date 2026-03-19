/*
 * Software implementation of se.h functions.
 * PKE registers backed by memory; Montgomery multiplication in software.
 */

#include "se.h"
#include <string.h>

/* 48 PKE registers, each 64 bytes = 16 uint32_t words */
#define PKE_NUM_REGS  48
#define PKE_REG_WORDS 16

static uint32_t pke_regs[PKE_NUM_REGS][PKE_REG_WORDS];

/* ---- No-op functions ---- */

int rsaPKESetLen(int len) {
  (void)len;
  return 0;
}

int eccCalMc64(uint32_t *mc, const uint32_t *p) {
  (void)mc;
  (void)p;
  return 0;
}

int rsaPKEWriteMc(const uint32_t *buf) {
  (void)buf;
  return 0;
}

/* ---- PKE register read/write ---- */

int eccPKEWriteBuf(int idx, const uint32_t *buf, int len) {
  if (idx < 0 || idx >= PKE_NUM_REGS || len < 0 || len > PKE_REG_WORDS)
    return -1;
  memset(pke_regs[idx], 0, sizeof(pke_regs[idx]));
  memcpy(pke_regs[idx], buf, (size_t)len * sizeof(uint32_t));
  return 0;
}

int eccPKEReadBuf(uint32_t *buf, int idx, int len) {
  if (idx < 0 || idx >= PKE_NUM_REGS || len < 0 || len > PKE_REG_WORDS)
    return -1;
  memcpy(buf, pke_regs[idx], (size_t)len * sizeof(uint32_t));
  return 0;
}

/* ---- 32-bit Montgomery multiplication ---- */

/* Cached Montgomery constant: -mod^{-1} mod 2^32 */
static uint32_t cached_mod = 0;
static uint32_t cached_neg_inv = 0;

/* Compute n^{-1} mod 2^32 using Newton's method (n must be odd) */
static uint32_t mod_inv32(uint32_t n) {
  uint32_t x = n;
  for (int i = 0; i < 5; i++)
    x *= 2 - n * x;
  return x;
}

static void ensure_mont_const(uint32_t mod) {
  if (mod != cached_mod) {
    cached_mod = mod;
    cached_neg_inv = -(mod_inv32(mod));
  }
}

int sm2MonMul(int r, int a, int b) {
  uint32_t mod = pke_regs[0][0];
  ensure_mont_const(mod);
  uint64_t t = (uint64_t)pke_regs[a][0] * pke_regs[b][0];
  uint32_t m = (uint32_t)t * cached_neg_inv;
  uint64_t u = t + (uint64_t)m * mod;
  uint32_t result = (uint32_t)(u >> 32);
  if (result >= mod) result -= mod;
  pke_regs[r][0] = result;
  return 0;
}

/* ---- Modular add/sub ---- */

int sm2MonAdd(int r, int a, int b) {
  (void)r;
  (void)a;
  (void)b;
  return 0;
}

int sm2MonSub(int r, int a, int b) {
  (void)r;
  (void)a;
  (void)b;
  return 0;
}

/* ---- Stubs ---- */

void self_test(void) {}
void se_start(void) {}
void se_stop(void) {}

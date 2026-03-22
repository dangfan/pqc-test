/*
 * ML-DSA-65 (FIPS 204)
 *
 * Design constraints
 * ------------------
 * Stack budget: total stack usage of any call path must stay under 5 KB.
 * Use PKE registers (48 × 64 B = 3 KB) to spill / cache polynomials
 */

#include "ml-dsa-65.h"
#include "se.h"
#include "sha3.h"
#include <string.h>
#include <stdio.h>
#include <stdint.h>
extern uint32_t uwTick;

#define Q       MLDSA_Q
#define N       MLDSA_N
#define K       MLDSA_K
#define L       MLDSA_L
#define ETA     MLDSA_ETA
#define TAU     MLDSA_TAU
#define BETA_B  MLDSA_BETA
#define GAMMA1  MLDSA_GAMMA1
#define GAMMA2  MLDSA_GAMMA2
#define OMEGA   MLDSA_OMEGA
#define D_BITS  MLDSA_D

#define QINV     58728449u        /* Q^{-1} mod 2^32 */
#define INTT_F   41978            /* R^2 * 256^{-1} mod Q  (compensates R^{-1} from pointwise) */

/* Kronecker substitution parameters */
#define KRON_G       16    /* coefficients per group */
#define KRON_T       16    /* number of groups (N / KRON_G) */
#define KRON_L       50    /* bits per packed coefficient */
#define KRON_PACK_WORDS  25  /* ceil(KRON_G * KRON_L / 32) = 800/32 */
#define KRON_PROD_WORDS  50  /* ceil((2*KRON_G - 1) * KRON_L / 32) ≈ 1550/32 */
#define KRON_PKE_LEN     54  /* pkeLen in 4-byte words for big-int multiply (>= KRON_PROD_WORDS + margin) */

typedef int32_t poly[N];

/* PKE register helpers: 16 int32 coefficients per 64-byte register */
#define PKE_COEFFS_PER_REG 16
#define PKE_SLOT0  4   /* regs  4-19: poly slot 0 */
#define PKE_SLOT1 20   /* regs 20-35: poly slot 1 */

/* 3-slot layout using all 48 regs (for sign path — no Kronecker needed) */
#define PKE_ACC0   0   /* regs  0-15: accumulator for row batch slot 0 */
#define PKE_ACC1  16   /* regs 16-31: accumulator for row batch slot 1 */
#define PKE_ACC2  32   /* regs 32-47: accumulator for row batch slot 2 */

/* NTT zetas in Montgomery domain (zeta_i * 2^32 mod Q).
 * Generated: zetas[i] = R * pow(1753, bitrev8(i), Q) centered. */
static const int32_t zetas[N] = {
 -4186625,   25847,-2608894, -518909,  237124, -777960, -876248,  466468,
  1826347, 2353451, -359251,-2091905, 3119733,-2884855, 3111497, 2680103,
  2725464, 1024112,-1079900, 3585928, -549488,-1119584, 2619752,-2108549,
 -2118186,-3859737,-1399561,-3277672, 1757237,  -19422, 4010497,  280005,
  2706023,   95776, 3077325, 3530437,-1661693,-3592148,-2537516, 3915439,
 -3861115,-3043716, 3574422,-2867647, 3539968, -300467, 2348700, -539299,
 -1699267,-1643818, 3505694,-3821735, 3507263,-2140649,-1600420, 3699596,
   811944,  531354,  954230, 3881043, 3900724,-2556880, 2071892,-2797779,
 -3930395,-1528703,-3677745,-3041255,-1452451, 3475950, 2176455,-1585221,
 -1257611, 1939314,-4083598,-1000202,-3190144,-3157330,-3632928,  126922,
  3412210, -983419, 2147896, 2715295,-2967645,-3693493, -411027,-2477047,
  -671102,-1228525,  -22981,-1308169, -381987, 1349076, 1852771,-1430430,
 -3343383,  264944,  508951, 3097992,   44288,-1100098,  904516, 3958618,
 -3724342,   -8578, 1653064,-3249728, 2389356, -210977,  759969,-1316856,
   189548,-3553272, 3159746,-1851402,-2409325, -177440, 1315589, 1341330,
  1285669,-1584928, -812732,-1439742,-3019102,-3881060,-3628969, 3839961,
  2091667, 3407706, 2316500, 3817976,-3342478, 2244091,-2446433,-3562462,
   266997, 2434439,-1235728, 3513181,-3520352,-3759364,-1197226,-3193378,
   900702, 1859098,  909542,  819034,  495491,-1613174,  -43260, -522500,
  -655327,-3122442, 2031748, 3207046,-3556995, -525098, -768622,-3595838,
   342297,  286988,-2437823, 4108315, 3437287,-3342277, 1735879,  203044,
  2842341, 2691481,-2590150, 1265009, 4055324, 1247620, 2486353, 1595974,
 -3767016, 1250494, 2635921,-3548272,-2994039, 1869119, 1903435,-1050970,
 -1333058, 1237275,-3318210,-1430225, -451100, 1312455, 3306115,-1962642,
 -1279661, 1917081,-2546312,-1374803, 1500165,  777191, 2235880, 3406031,
  -542412,-2831860,-1671176,-1846953,-2584293,-3724270,  594136,-3776993,
 -2013608, 2432395, 2454455, -164721, 1957272, 3369112,  185531,-1207385,
 -3183426,  162844, 1616392, 3014001,  810149, 1652634,-3694233,-1799107,
 -3038916, 3523897, 3866901,  269760, 2213111, -975884, 1717735,  472078,
  -426683, 1723600,-1803090, 1910376,-1667432,-1104333, -260646,-3833893,
 -2939036,-2235985, -420899,-2286327,  183443, -976891, 1612842,-3545687,
  -554416, 3919660,  -48306,-1362209, 3937738, 1400424, -846154, 1976782
};

/* ------------------------------------------------------------------ */
/*  Modular arithmetic helpers                                        */
/* ------------------------------------------------------------------ */

/* Dilithium-style fast reduction: avoids software division.
 * reduce32: map a to range (-Q/2, Q/2) via shift-multiply.
 * caddq:    conditional add Q if negative.
 * freeze:   full reduction to [0, Q-1]. */
static inline int32_t reduce32(int32_t a) {
  int32_t t = (a + (1 << 22)) >> 23;
  return a - t * Q;
}

static inline int32_t caddq(int32_t a) {
  return a + ((a >> 31) & Q);
}

static int32_t freeze(int32_t a) {
  return caddq(reduce32(a));
}

/* ------------------------------------------------------------------ */
/*  PKE register read/write for polynomial storage                    */
/*  Store coefficients as native uint32_t words (16 per register).    */
/* ------------------------------------------------------------------ */

static void pke_store_poly(int base_reg, const int32_t p[N]) {
  for (int r = 0; r < N / PKE_COEFFS_PER_REG; r++)
    eccPKEWriteBuf(base_reg + r, (const uint32_t *)&p[r * PKE_COEFFS_PER_REG],
                   PKE_COEFFS_PER_REG);
}

static void pke_load_poly(int32_t p[N], int base_reg) {
  for (int r = 0; r < N / PKE_COEFFS_PER_REG; r++)
    eccPKEReadBuf((uint32_t *)&p[r * PKE_COEFFS_PER_REG], base_reg + r,
                  PKE_COEFFS_PER_REG);
}

/* Software Montgomery multiply: a * b * 2^{-32} mod Q.
 * Standard Dilithium Montgomery reduction. */
static int32_t montgomery_reduce(int64_t a) {
  int32_t t;
  t = (int64_t)(int32_t)a * QINV;
  t = (a - (int64_t)t * Q) >> 32;
  return t;
}

static int32_t mon_mul(int32_t a, int32_t b) {
  return montgomery_reduce((int64_t)a * b);
}

/* ------------------------------------------------------------------ */
/*  NTT / INTT (Cooley-Tukey / Gentleman-Sande, software Montgomery)  */
/* ------------------------------------------------------------------ */

static void ntt(int32_t a[N]) {
  unsigned int len, start, j, k;
  int32_t t;

  k = 0;
  for (len = 128; len >= 1; len >>= 1) {
    for (start = 0; start < N; start = j + len) {
      int32_t zeta = zetas[++k];
      for (j = start; j < start + len; ++j) {
        t = mon_mul(zeta, a[j + len]);
        a[j + len] = a[j] - t;
        a[j]       = a[j] + t;
      }
    }
  }
}

static void invntt(int32_t a[N]) {
  unsigned int start, len, j, k;
  int32_t t;

  /* Reduce inputs to (-Q, Q) to prevent int32 overflow in butterfly.
     The upper position a[j] = t + a[j+len] doubles each stage;
     after 8 stages the worst case is 256*Q ≈ 2.15B < INT32_MAX. */
  for (j = 0; j < N; j++)
    a[j] = reduce32(a[j]);

  k = 256;
  for (len = 1; len <= 128; len <<= 1) {
    for (start = 0; start < N; start = j + len) {
      int32_t zeta = -zetas[--k];
      for (j = start; j < start + len; ++j) {
        t = a[j];
        a[j]       = t + a[j + len];
        a[j + len] = t - a[j + len];
        a[j + len] = mon_mul(zeta, a[j + len]);
      }
    }
  }
  for (j = 0; j < N; ++j)
    a[j] = mon_mul(INTT_F, a[j]);
}

/* c += a * b (NTT domain, Montgomery, accumulate) */
static void poly_pointwise_acc(int32_t c[N], const int32_t a[N],
                               const int32_t b[N]) {
  for (int i = 0; i < N; i++)
    c[i] += mon_mul(a[i], b[i]);
}

/* -------------------------------------------------------------------- */
/*  Kronecker substitution polynomial multiplication                    */
/*                                                                      */
/*  Replaces NTT + pointwise + INTT with direct polynomial multiply     */
/*  in Z_q[X]/(X^N+1) using Kronecker packing and hardware big-int MMM. */
/*                                                                      */
/*  Strategy:                                                           */
/*    - Split each poly into KRON_T=16 groups of KRON_G=16 coefficients */
/*    - For each pair of groups (i,j), pack coefficients at KRON_L=50   */
/*      bits each into ~800-bit integers, do ONE hardware multiply,     */
/*      unpack the 31-coefficient sub-product.                          */
/*    - Accumulate sub-products into result, reducing mod X^N+1.        */
/*    - Total: 256 HW multiplies (vs ~2304 in NTT approach).            */
/* -------------------------------------------------------------------- */

/* Pack KRON_G non-negative coefficients into a multi-word integer.
 * Each coefficient occupies KRON_L bits at a KRON_L-bit stride.
 * Result is in host word order (word[0] = least significant). */
static void kron_pack(uint32_t packed[KRON_PACK_WORDS],
                      const int32_t *coeffs, int ncoeffs) {
  memset(packed, 0, KRON_PACK_WORDS * 4);
  for (int i = 0; i < ncoeffs; i++) {
    uint64_t val = (uint64_t)(uint32_t)coeffs[i];  /* non-negative */
    int bit_pos = i * KRON_L;
    int word = bit_pos / 32;
    int shift = bit_pos % 32;
    packed[word] |= (uint32_t)(val << shift);
    if (shift + KRON_L > 32 && word + 1 < KRON_PACK_WORDS)
      packed[word + 1] |= (uint32_t)(val >> (32 - shift));
    if (shift + KRON_L > 64 && word + 2 < KRON_PACK_WORDS)
      packed[word + 2] |= (uint32_t)(val >> (64 - shift));
  }
}

/* Unpack a single coefficient at index k from a product big integer. */
static int64_t kron_unpack_one(int k, const uint32_t prod[KRON_PROD_WORDS]) {
  uint64_t mask = ((uint64_t)1 << KRON_L) - 1;
  int bit_pos = k * KRON_L;
  int word = bit_pos / 32;
  int shift = bit_pos % 32;
  uint64_t val = (uint64_t)prod[word] >> shift;
  if (shift + KRON_L > 32 && word + 1 < KRON_PROD_WORDS)
    val |= (uint64_t)prod[word + 1] << (32 - shift);
  if (shift + KRON_L > 64 && word + 2 < KRON_PROD_WORDS)
    val |= (uint64_t)prod[word + 2] << (64 - shift);
  return (int64_t)(val & mask);
}

/* Big-integer MMM setup for Kronecker multiplication.
 * Uses modulus N = 2^M - 1 (Mersenne-like) so R ≡ 1 (mod N),
 * making MonMul(a,b) = a*b mod N = a*b when a*b < N. */
static void kron_mmm_setup(void) {
  rsaPKESetLen(KRON_PKE_LEN);

  /* Set modulus N = 2^(KRON_PKE_LEN*32) - 1 (all ones, odd).
   * With N = 2^M - 1 where M = KRON_PKE_LEN*32:
   *   R = 2^M ≡ 1 (mod N), so R^{-1} ≡ 1 (mod N).
   *   MonMul(a, b) = a * b * R^{-1} mod N = a * b mod N.
   * As long as a*b < N, we get exact a*b. No R compensation needed. */
  uint32_t mod[KRON_PKE_LEN];
  memset(mod, 0xFF, KRON_PKE_LEN * 4);
  eccPKEWriteBuf(0, mod, KRON_PKE_LEN);

  uint32_t mc[2];
  eccCalMc64(mc, mod);
  rsaPKEWriteMc(mc);
}

/* Convert a multi-word integer between host representation and PKE
 * register representation, in-place.
 *
 * Host representation: uint32_t words in native byte order, word[0] is
 * the least significant (least-significant-word-first).
 *
 * PKE representation: bytes in big-endian order, i.e. the most
 * significant byte at the lowest address (most-significant-word-first,
 * each word stored big-endian).
 *
 * The transform is self-inverse: applying it twice restores the original. */
static void host_to_pke(uint32_t *buf, int nwords) {
  int i = 0, j = nwords - 1;
  for (; i < j; i++, j--) {
    uint32_t a = buf[i], b = buf[j];
#if __BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__
    buf[i] = __builtin_bswap32(b);
    buf[j] = __builtin_bswap32(a);
#else
    buf[i] = b;
    buf[j] = a;
#endif
  }
#if __BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__
  if (i == j)
    buf[i] = __builtin_bswap32(buf[i]);
#endif
}
/* pke_to_host is the same operation (self-inverse). */
#define pke_to_host(buf, nwords) host_to_pke((buf), (nwords))

/* Multiply two packed big integers using hardware MMM.
 * Since N = 2^M - 1 and R = 2^M ≡ 1 (mod N), MonMul(a,b) = a*b mod N.
 * As long as a*b < N, we get exact a*b.
 *
 * PKE registers store big-endian multi-precision integers.
 * Packed arrays use host uint32_t words in LSW-first order.
 * host_to_pke / pke_to_host convert between the two. */
/* Register layout for Kronecker big-int multiply (pkeLen=54, 216 bytes each):
 * Each operand spans ceil(216/64) = 4 register slots.
 * Reg  0- 3: modulus (set by kron_mmm_setup)
 * Reg 36-39: operand a    (avoids PKE_SLOT0=4-19, PKE_SLOT1=20-35)
 * Reg 40-43: operand b
 * Reg 44-47: result */
#define KRON_REG_A  36
#define KRON_REG_B  40
#define KRON_REG_R  44

/* Multiply challenge (in PKE_SLOT0) by polynomial b in Z_q[X]/(X^N+1).
 * Challenge c is read from PKE registers — no stack poly needed for it.
 * Accumulates into int32_t result with per-step mod Q (no int64_t acc[N]).
 * b must have coefficients in [0, Q). Result in [0, Q).
 *
 * Kronecker big-int multiply is inlined to share the tmp buffer for both
 * operand writes and product readback, eliminating a separate prod[] array
 * and the kron_bigint_mul call frame.  The 'a' operand is written to PKE
 * once per outer (gi) iteration instead of every inner (gj) iteration. */
static void poly_mul_challenge(int32_t result[N], const int32_t b[N]) {
  int32_t a_chunk[KRON_G];               /* 64 bytes: challenge group */
  uint32_t a_packed[KRON_PACK_WORDS];    /* 100 bytes */
  uint32_t b_packed[KRON_PACK_WORDS];    /* 100 bytes */
  uint32_t tmp[KRON_PKE_LEN];            /* 216 bytes: shared PKE I/O + product */

  memset(result, 0, N * sizeof(int32_t));
  kron_mmm_setup();

  /* Inlined Kronecker big-int multiply (cf. kron_bigint_mul).
   *
   * Compared to calling kron_bigint_mul per (gi, gj) pair:
   *   1. The shared tmp[KRON_PKE_LEN] serves triple duty — write operand a
   *      to PKE, write operand b to PKE, and hold the product readback —
   *      eliminating a separate prod[KRON_PROD_WORDS] (200 B) on the stack.
   *   2. The kron_bigint_mul call frame (~288 B) is eliminated entirely,
   *      flattening the deepest call chain of the sign path.
   *   3. Operand 'a' (the challenge group) is written to KRON_REG_A once
   *      per outer gi iteration; the original code rewrote it every (gi,gj)
   *      pair — saving 16 × (KRON_T − 1) redundant PKE writes per call. */
  for (int gi = 0; gi < KRON_T; gi++) {
    eccPKEReadBuf((uint32_t *)a_chunk, PKE_SLOT0 + gi, KRON_G);
    kron_pack(a_packed, a_chunk, KRON_G);

    /* Pack challenge group into big-endian PKE format and write once */
    memset(tmp, 0, KRON_PKE_LEN * 4);
    memcpy(tmp, a_packed, KRON_PACK_WORDS * 4);
    host_to_pke(tmp, KRON_PKE_LEN);
    eccPKEWriteBuf(KRON_REG_A, tmp, KRON_PKE_LEN);

    for (int gj = 0; gj < KRON_T; gj++) {
      kron_pack(b_packed, &b[gj * KRON_G], KRON_G);

      /* Pack b group and write to PKE */
      memset(tmp, 0, KRON_PKE_LEN * 4);
      memcpy(tmp, b_packed, KRON_PACK_WORDS * 4);
      host_to_pke(tmp, KRON_PKE_LEN);
      eccPKEWriteBuf(KRON_REG_B, tmp, KRON_PKE_LEN);

      /* Hardware big-int multiply: KRON_REG_R = a * b mod (2^M − 1) */
      sm2MonMul(KRON_REG_R, KRON_REG_A, KRON_REG_B);

      /* Reuse tmp for product readback (host word order) */
      eccPKEReadBuf(tmp, KRON_REG_R, KRON_PKE_LEN);
      pke_to_host(tmp, KRON_PKE_LEN);

      int base = (gi + gj) * KRON_G;
      for (int k = 0; k < 2 * KRON_G - 1; k++) {
        int64_t coeff = kron_unpack_one(k, tmp);
        int idx = base + k;
        int64_t v;
        if (idx < N) {
          v = (int64_t)result[idx] + coeff;
        } else {
          idx -= N;
          v = (int64_t)result[idx] - coeff;
        }
        v %= Q;
        if (v < 0) v += Q;
        result[idx] = (int32_t)v;
      }
    }
  }
}

static void poly_add(int32_t c[N], const int32_t a[N], const int32_t b[N]) {
  for (int i = 0; i < N; i++) c[i] = a[i] + b[i];
}

static void poly_sub(int32_t c[N], const int32_t a[N], const int32_t b[N]) {
  for (int i = 0; i < N; i++) c[i] = a[i] - b[i];
}

static void poly_reduce(int32_t a[N]) {
  for (int i = 0; i < N; i++) a[i] = freeze(a[i]);
}

static void poly_caddq(int32_t a[N]) {
  for (int i = 0; i < N; i++) a[i] = caddq(a[i]);
}

/* ------------------------------------------------------------------ */
/*  Sampling: RejNTTPoly (ExpandA), ExpandMask, RejBoundedPoly,       */
/*            SampleInBall                                            */
/* ------------------------------------------------------------------ */

/* Rejection-sample one NTT-domain polynomial from SHAKE128 (ExpandA).
 * FIPS 204, Algorithm 30 (RejNTTPoly). */
static void poly_rej_ntt(int32_t a[N], const uint8_t rho[32],
                         uint8_t row, uint8_t col) {
  SHA3_CTX_T ctx;
  /* Squeeze full SHAKE128 blocks (168 bytes = 56 samples of 3 bytes) */
  uint8_t buf[168];
  int ctr = 0;
  int bufpos = 168; /* start exhausted → trigger first squeeze */

  shake128_init(&ctx);
  {
    uint8_t hdr[2] = {col, row};
    shake_update(&ctx, rho, 32);
    shake_update(&ctx, hdr, 2);
  }
  shake_finalize(&ctx);

  while (ctr < N) {
    if (bufpos >= 168) {
      shake_squeeze(&ctx, buf, 168);
      bufpos = 0;
    }
    uint32_t t = ((uint32_t)buf[bufpos]) | ((uint32_t)buf[bufpos + 1] << 8) |
                 ((uint32_t)(buf[bufpos + 2] & 0x7F) << 16);
    bufpos += 3;
    if (t < (uint32_t)Q)
      a[ctr++] = (int32_t)t;
  }
}

/* Expand masking vector component y_idx from rho'.
 * FIPS 204, Algorithm 36 (ExpandMask). gamma1 = 2^19, 20-bit packing. */
static void poly_expand_mask(int32_t a[N], const uint8_t rho_prime[64],
                             uint16_t idx) {
  SHA3_CTX_T ctx;
  /* N/2 pairs × 5 bytes = 640 bytes total, no rejection needed */
  uint8_t buf[N / 2 * 5];
  uint8_t hdr[2];

  shake256_init(&ctx);
  shake_update(&ctx, rho_prime, 64);
  hdr[0] = (uint8_t)(idx & 0xFF);
  hdr[1] = (uint8_t)(idx >> 8);
  shake_update(&ctx, hdr, 2);
  shake_finalize(&ctx);

  shake_squeeze(&ctx, buf, sizeof(buf));

  /* 20-bit packed: 2 coefficients per 5 bytes */
  for (int i = 0; i < N / 2; i++) {
    int off = i * 5;
    uint32_t t0 = ((uint32_t)buf[off]) | ((uint32_t)buf[off + 1] << 8) |
                  ((uint32_t)(buf[off + 2] & 0x0F) << 16);
    uint32_t t1 = ((uint32_t)buf[off + 2] >> 4) | ((uint32_t)buf[off + 3] << 4) |
                  ((uint32_t)buf[off + 4] << 12);
    a[2 * i]     = (int32_t)(GAMMA1 - t0);
    a[2 * i + 1] = (int32_t)(GAMMA1 - t1);
  }
}

/* Rejection-sample a bounded polynomial with coefficients in [-eta, eta].
 * FIPS 204, Algorithm 31 (RejBoundedPoly). eta=4: sample nibbles. */
static void poly_rej_bounded(int32_t a[N], const uint8_t seed[64],
                             uint16_t nonce) {
  SHA3_CTX_T ctx;
  /* Squeeze full SHAKE256 blocks (136 bytes) at a time */
  uint8_t buf[136];
  uint8_t hdr[2];
  int ctr = 0;
  int bufpos = 136; /* start exhausted → trigger first squeeze */

  shake256_init(&ctx);
  shake_update(&ctx, seed, 64);
  hdr[0] = (uint8_t)(nonce & 0xFF);
  hdr[1] = (uint8_t)(nonce >> 8);
  shake_update(&ctx, hdr, 2);
  shake_finalize(&ctx);

  while (ctr < N) {
    if (bufpos >= 136) {
      shake_squeeze(&ctx, buf, 136);
      bufpos = 0;
    }
    uint8_t z0 = buf[bufpos] & 0x0F;
    uint8_t z1 = buf[bufpos] >> 4;
    bufpos++;
    if (z0 <= 2 * ETA) {
      a[ctr++] = (int32_t)ETA - (int32_t)z0;
      if (ctr >= N) break;
    }
    if (z1 <= 2 * ETA) {
      a[ctr++] = (int32_t)ETA - (int32_t)z1;
    }
  }
}

/* SampleInBall: produce challenge polynomial c with TAU +/-1 coefficients.
 * FIPS 204, Algorithm 32. */
static void poly_challenge(int32_t c[N], const uint8_t c_tilde[48]) {
  SHA3_CTX_T ctx;
  /* Squeeze full SHAKE256 blocks (136 bytes) at a time */
  uint8_t buf[136];
  uint64_t signs;
  int bufpos;

  memset(c, 0, N * sizeof(int32_t));

  shake256_init(&ctx);
  shake_update(&ctx, c_tilde, MLDSA_C_TILDE_BYTES);
  shake_finalize(&ctx);
  shake_squeeze(&ctx, buf, 136);
  signs = 0;
  for (int i = 0; i < 8; i++)
    signs |= (uint64_t)buf[i] << (8 * i);
  bufpos = 8;

  for (unsigned int i = N - TAU; i < N; i++) {
    uint8_t b;
    do {
      if (bufpos >= 136) {
        shake_squeeze(&ctx, buf, 136);
        bufpos = 0;
      }
      b = buf[bufpos++];
    } while ((unsigned int)b > i);
    unsigned int j = (unsigned int)b;
    c[i] = c[j];
    c[j] = (signs & 1) ? -1 : 1;
    signs >>= 1;
  }
}

/* ------------------------------------------------------------------ */
/*  Sparse challenge multiplication                                    */
/*                                                                     */
/*  The challenge c has only TAU=49 non-zero coefficients, each ±1.    */
/*  Instead of dense Kronecker multiply (256 HW sm2MonMul), we do      */
/*  TAU passes of N add/sub — pure software, ~50× faster.             */
/* ------------------------------------------------------------------ */

typedef struct {
  uint8_t pos[TAU];   /* indices of non-zero coefficients */
  int8_t  sign[TAU];  /* +1 or -1 */
} challenge_sparse_t;

/* Extract sparse representation from challenge polynomial (values ±1). */
static void challenge_to_sparse(challenge_sparse_t *ch, const int32_t c[N]) {
  int k = 0;
  for (int i = 0; i < N && k < TAU; i++) {
    if (c[i] == 1) {
      ch->pos[k] = (uint8_t)i;
      ch->sign[k] = 1;
      k++;
    } else if (c[i] == -1) {
      ch->pos[k] = (uint8_t)i;
      ch->sign[k] = -1;
      k++;
    }
  }
}

/* Centered sparse multiply: b[] is centered (e.g. [-ETA,ETA] or
 * [-(2^12-1), 2^12]).  Output is centered with NO mod-Q reduction.
 *
 * Range analysis (worst case, all TAU taps align on one coefficient):
 *   eta-bounded b ∈ [-4, 4]:       |result[i]| ≤ TAU*ETA   =  196
 *   t0-bounded  b ∈ [-4095, 4096]: |result[i]| ≤ TAU*4096  = 200704
 * Both fit comfortably in int32_t — no intermediate reduction needed. */
static void poly_mul_sparse_centered(int32_t result[N], const int32_t b[N],
                                     const challenge_sparse_t *ch) {
  memset(result, 0, N * sizeof(int32_t));
  for (int k = 0; k < TAU; k++) {
    int j = ch->pos[k];
    int s = ch->sign[k];
    if (s > 0) {
      for (int i = 0; i < N - j; i++)
        result[i + j] += b[i];
      for (int i = N - j; i < N; i++)
        result[i + j - N] -= b[i];
    } else {
      for (int i = 0; i < N - j; i++)
        result[i + j] -= b[i];
      for (int i = N - j; i < N; i++)
        result[i + j - N] += b[i];
    }
  }
  /* No reduction — output stays centered. */
}

/* ------------------------------------------------------------------ */
/*  Packing / Unpacking                                               */
/* ------------------------------------------------------------------ */

/* Unpack eta-bounded poly from secret key (eta=4, 4 bits per coeff). */
static void unpack_eta(int32_t a[N], const uint8_t *buf) {
  for (int i = 0; i < N / 2; i++) {
    uint8_t b = buf[i];
    a[2 * i]     = (int32_t)ETA - (int32_t)(b & 0x0F);
    a[2 * i + 1] = (int32_t)ETA - (int32_t)(b >> 4);
  }
}

/* Unpack t0 polynomial (13 bits per coefficient, unsigned [0, 2^13)). */
static void unpack_t0(int32_t a[N], const uint8_t *buf) {
  for (int i = 0; i < N / 8; i++) {
    const uint8_t *p = buf + i * 13;
    a[8*i+0] = (1 << 12) - (int32_t)(((uint32_t)p[0] | ((uint32_t)p[1] << 8)) & 0x1FFF);
    a[8*i+1] = (1 << 12) - (int32_t)((((uint32_t)p[1] >> 5) | ((uint32_t)p[2] << 3) | ((uint32_t)p[3] << 11)) & 0x1FFF);
    a[8*i+2] = (1 << 12) - (int32_t)((((uint32_t)p[3] >> 2) | ((uint32_t)p[4] << 6)) & 0x1FFF);
    a[8*i+3] = (1 << 12) - (int32_t)((((uint32_t)p[4] >> 7) | ((uint32_t)p[5] << 1) | ((uint32_t)p[6] << 9)) & 0x1FFF);
    a[8*i+4] = (1 << 12) - (int32_t)((((uint32_t)p[6] >> 4) | ((uint32_t)p[7] << 4) | ((uint32_t)p[8] << 12)) & 0x1FFF);
    a[8*i+5] = (1 << 12) - (int32_t)((((uint32_t)p[8] >> 1) | ((uint32_t)p[9] << 7)) & 0x1FFF);
    a[8*i+6] = (1 << 12) - (int32_t)((((uint32_t)p[9] >> 6) | ((uint32_t)p[10] << 2) | ((uint32_t)p[11] << 10)) & 0x1FFF);
    a[8*i+7] = (1 << 12) - (int32_t)((((uint32_t)p[11] >> 3) | ((uint32_t)p[12] << 5)) & 0x1FFF);
  }
}

/* Pack z polynomial (gamma1=2^19, 20 bits per coeff). */
static void pack_z(uint8_t *buf, const int32_t a[N]) {
  for (int i = 0; i < N / 2; i++) {
    uint32_t t0 = (uint32_t)(GAMMA1 - a[2 * i]);
    uint32_t t1 = (uint32_t)(GAMMA1 - a[2 * i + 1]);
    buf[5*i+0] = (uint8_t)(t0);
    buf[5*i+1] = (uint8_t)(t0 >> 8);
    buf[5*i+2] = (uint8_t)((t0 >> 16) | (t1 << 4));
    buf[5*i+3] = (uint8_t)(t1 >> 4);
    buf[5*i+4] = (uint8_t)(t1 >> 12);
  }
}

/* Pack t1 polynomial (10 bits per coeff). */
static void pack_t1(uint8_t *buf, const int32_t a[N]) {
  for (int i = 0; i < N / 4; i++) {
    buf[5*i+0] = (uint8_t)(a[4*i+0]);
    buf[5*i+1] = (uint8_t)((a[4*i+0] >> 8) | (a[4*i+1] << 2));
    buf[5*i+2] = (uint8_t)((a[4*i+1] >> 6) | (a[4*i+2] << 4));
    buf[5*i+3] = (uint8_t)((a[4*i+2] >> 4) | (a[4*i+3] << 6));
    buf[5*i+4] = (uint8_t)(a[4*i+3] >> 2);
  }
}

/* Pack eta-bounded poly (eta=4, 4 bits per coeff). */
static void pack_eta(uint8_t *buf, const int32_t a[N]) {
  for (int i = 0; i < N / 2; i++) {
    uint8_t v0 = (uint8_t)(ETA - a[2*i]);
    uint8_t v1 = (uint8_t)(ETA - a[2*i+1]);
    buf[i] = (v0 & 0x0F) | (v1 << 4);
  }
}

/* Pack t0 polynomial (13 bits per coeff). */
static void pack_t0(uint8_t *buf, const int32_t a[N]) {
  for (int i = 0; i < N / 8; i++) {
    uint32_t t[8];
    for (int k = 0; k < 8; k++)
      t[k] = (uint32_t)((1 << 12) - a[8*i+k]);
    buf[13*i+ 0] = (uint8_t)(t[0]);
    buf[13*i+ 1] = (uint8_t)((t[0] >> 8) | (t[1] << 5));
    buf[13*i+ 2] = (uint8_t)(t[1] >> 3);
    buf[13*i+ 3] = (uint8_t)((t[1] >> 11) | (t[2] << 2));
    buf[13*i+ 4] = (uint8_t)((t[2] >> 6) | (t[3] << 7));
    buf[13*i+ 5] = (uint8_t)(t[3] >> 1);
    buf[13*i+ 6] = (uint8_t)((t[3] >> 9) | (t[4] << 4));
    buf[13*i+ 7] = (uint8_t)(t[4] >> 4);
    buf[13*i+ 8] = (uint8_t)((t[4] >> 12) | (t[5] << 1));
    buf[13*i+ 9] = (uint8_t)((t[5] >> 7) | (t[6] << 6));
    buf[13*i+10] = (uint8_t)(t[6] >> 2);
    buf[13*i+11] = (uint8_t)((t[6] >> 10) | (t[7] << 3));
    buf[13*i+12] = (uint8_t)(t[7] >> 5);
  }
}

/* ----------------------------------------------------------------- */
/*  Decompose, HighBits, LowBits, MakeHint, Power2Round              */
/* ----------------------------------------------------------------- */

/* Decompose: r = r1 * 2*GAMMA2 + r0, with r0 in (-GAMMA2, GAMMA2].
 * For gamma2 = (q-1)/32, r1 in [0, 15]. */
static int32_t decompose(int32_t *r0, int32_t a) {
  int32_t a1 = (a + 127) >> 7;
  a1 = (a1 * 1025 + (1 << 21)) >> 22;
  a1 &= 15;
  *r0 = a - a1 * 2 * GAMMA2;
  *r0 -= (((Q - 1) / 2 - *r0) >> 31) & Q;
  return a1;
}

static int32_t high_bits(int32_t a) {
  int32_t r0;
  return decompose(&r0, a);
}

/* Encode w1 = HighBits(w) and feed directly into a SHAKE context.
 * Processes 32 coefficients (16 bytes) at a time, avoiding a full
 * POLYW1_PACKEDBYTES (128-byte) intermediate buffer on the stack. */
static void __attribute__((noinline)) w1_encode_update(SHA3_CTX_T *ctx,
                                                       const int32_t w[N]) {
  uint8_t chunk[16];
  for (int i = 0; i < N; i += 32) {
    for (int j = 0; j < 16; j++)
      chunk[j] = (uint8_t)((uint32_t)high_bits(w[i + 2*j])
                          | ((uint32_t)high_bits(w[i + 2*j + 1]) << 4));
    shake_update(ctx, chunk, 16);
  }
}

/* Pure-threshold hint from already-decomposed r and centered ct0.
 * Equivalent to: high_bits(r) != high_bits(freeze(r + ct0))
 * where r1 = decompose(&r0, r).
 *
 * Adding ct0 to r shifts the low part by ct0.  The high bits change
 * iff the shifted low part leaves the half-open interval (-GAMMA2, GAMMA2]:
 *   s = r0 + ct0
 *   hint = (s > GAMMA2) || (s < -GAMMA2) || (s == -GAMMA2 && r1 != 0)
 *
 * The r1 != 0 guard handles the wrap: for r1 == 0 the "negative"
 * representatives (r0 < 0, a near Q) stay in bin 0 when s == -GAMMA2,
 * while for r1 >= 1 the value crosses into the adjacent lower bin. */
static inline unsigned int make_hint_from_decomposed(int32_t r1,
                                                      int32_t r0,
                                                      int32_t ct0) {
  int32_t s = r0 + ct0;
  return (s > (int32_t)GAMMA2) | (s < -(int32_t)GAMMA2) |
         ((s == -(int32_t)GAMMA2) & (r1 != 0));
}

/* Power2Round: a = a1 * 2^D + a0 */
static int32_t power2round(int32_t *a0, int32_t a) {
  int32_t a1 = (a + (1 << (D_BITS - 1)) - 1) >> D_BITS;
  *a0 = a - (a1 << D_BITS);
  return a1;
}

/* ------------------------------------------------------------------ */
/*  Secret key layout accessors                                       */
/* ------------------------------------------------------------------ */

static const uint8_t *sk_rho(const uint8_t *sk)  { return sk; }
static const uint8_t *sk_K(const uint8_t *sk)    { return sk + 32; }
static const uint8_t *sk_tr(const uint8_t *sk)   { return sk + 64; }

static const uint8_t *sk_s1(const uint8_t *sk, int j) {
  return sk + 128 + j * MLDSA_POLYETA_PACKEDBYTES;
}
static const uint8_t *sk_s2(const uint8_t *sk, int i) {
  return sk + 128 + L * MLDSA_POLYETA_PACKEDBYTES
       + i * MLDSA_POLYETA_PACKEDBYTES;
}
static const uint8_t *sk_t0(const uint8_t *sk, int i) {
  return sk + 128 + (L + K) * MLDSA_POLYETA_PACKEDBYTES
       + i * MLDSA_POLYT0_PACKEDBYTES;
}

/* ------------------------------------------------------------------ */
/*  Compute w_hat[i] = sum_j A_hat[i][j] * y_hat[j]                   */
/*  One row at a time. y is re-expanded for each row.                 */
/*  Uses poly0 for y_hat_j, poly1 for A_hat_ij.                       */
/*  w_hat is accumulated in PKE_SLOT1 (no stack poly needed).         */
/*  Caller retrieves result via pke_load_poly(dst, PKE_SLOT1).        */
/* ------------------------------------------------------------------ */

/* Zero a 16-register PKE slot */
static void pke_zero_slot(int base_reg) {
  uint32_t zeros[PKE_COEFFS_PER_REG];
  memset(zeros, 0, sizeof(zeros));
  for (int r = 0; r < N / PKE_COEFFS_PER_REG; r++)
    eccPKEWriteBuf(base_reg + r, zeros, PKE_COEFFS_PER_REG);
}

/* Pointwise accumulate into any PKE slot: slot[k] += a[k] * b[k]. */
static void pointwise_acc_to(int base_reg,
                             const int32_t a[N], const int32_t b[N]) {
  int32_t w_chunk[PKE_COEFFS_PER_REG];
  for (int r = 0; r < N / PKE_COEFFS_PER_REG; r++) {
    eccPKEReadBuf((uint32_t *)w_chunk, base_reg + r, PKE_COEFFS_PER_REG);
    for (int k = 0; k < PKE_COEFFS_PER_REG; k++)
      w_chunk[k] += mon_mul(a[r * PKE_COEFFS_PER_REG + k],
                             b[r * PKE_COEFFS_PER_REG + k]);
    eccPKEWriteBuf(base_reg + r, (const uint32_t *)w_chunk,
                   PKE_COEFFS_PER_REG);
  }
}

/* Pointwise accumulate into PKE_SLOT1: w_hat += a * b (NTT domain).
 * Reads/writes w_hat from PKE_SLOT1 in 16-coeff chunks.
 * Only touches regs 0-3 (MMM) and 20-35 (PKE_SLOT1). */
static void pointwise_acc_pke(const int32_t a[N], const int32_t b[N]) {
  pointwise_acc_to(PKE_SLOT1, a, b);
}

/* Like pointwise_acc_pke, but reads b from PKE_SLOT0 instead of stack.
 * SLOT1 += a[k] * SLOT0[k].  Only needs one stack poly for a. */
static void pointwise_acc_pke_slot0(const int32_t a[N]) {
  int32_t w_chunk[PKE_COEFFS_PER_REG];
  int32_t b_chunk[PKE_COEFFS_PER_REG];
  for (int r = 0; r < N / PKE_COEFFS_PER_REG; r++) {
    eccPKEReadBuf((uint32_t *)w_chunk, PKE_SLOT1 + r, PKE_COEFFS_PER_REG);
    eccPKEReadBuf((uint32_t *)b_chunk, PKE_SLOT0 + r, PKE_COEFFS_PER_REG);
    for (int k = 0; k < PKE_COEFFS_PER_REG; k++)
      w_chunk[k] += mon_mul(a[r * PKE_COEFFS_PER_REG + k], b_chunk[k]);
    eccPKEWriteBuf(PKE_SLOT1 + r, (const uint32_t *)w_chunk,
                   PKE_COEFFS_PER_REG);
  }
}

/* Profiling accumulators — reset per iteration in ml_dsa_65_sign */
static uint32_t prof_expand_mask, prof_ntt_fwd, prof_rej_ntt, prof_pw_acc, prof_invntt;

static void compute_w_hat_row(const uint8_t *rho,
                              const uint8_t *rho_prime,
                              uint16_t kappa, int row,
                              int32_t poly0[N], int32_t poly1[N]) {
  /* Zero PKE_SLOT1 (w_hat accumulator) */
  pke_zero_slot(PKE_SLOT1);
  for (int j = 0; j < L; j++) {
    uint32_t _t;
    /* Expand y_j and NTT */
    _t = uwTick; poly_expand_mask(poly0, rho_prime, kappa + (uint16_t)j); prof_expand_mask += uwTick - _t;
    _t = uwTick; ntt(poly0); prof_ntt_fwd += uwTick - _t;
    /* Expand A[row][j] (already in NTT domain) */
    _t = uwTick; poly_rej_ntt(poly1, rho, (uint8_t)row, (uint8_t)j); prof_rej_ntt += uwTick - _t;
    /* w_hat += A_ij * y_hat_j (accumulated in PKE_SLOT1) */
    _t = uwTick; pointwise_acc_pke(poly1, poly0); prof_pw_acc += uwTick - _t;
  }
}

/* Compute w_hat for 3 rows simultaneously, expanding each y[j] only once.
 * Results in PKE_ACC0 (row0), PKE_ACC1 (row1), PKE_ACC2 (row2).
 * Reduces expand_mask+NTT calls from 3*L to L per batch of 3 rows. */
static void compute_w_hat_batch3(const uint8_t *rho,
                                 const uint8_t *rho_prime,
                                 uint16_t kappa, int row0,
                                 int32_t poly0[N], int32_t poly1[N]) {
  pke_zero_slot(PKE_ACC0);
  pke_zero_slot(PKE_ACC1);
  pke_zero_slot(PKE_ACC2);
  for (int j = 0; j < L; j++) {
    uint32_t _t;
    /* Expand y_j and NTT — done ONCE for all 3 rows */
    _t = uwTick; poly_expand_mask(poly0, rho_prime, kappa + (uint16_t)j); prof_expand_mask += uwTick - _t;
    _t = uwTick; ntt(poly0); prof_ntt_fwd += uwTick - _t;
    /* poly0 = y_hat[j], reused for 3 rows */
    /* Row 0 */
    _t = uwTick; poly_rej_ntt(poly1, rho, (uint8_t)(row0),     (uint8_t)j); prof_rej_ntt += uwTick - _t;
    _t = uwTick; pointwise_acc_to(PKE_ACC0, poly1, poly0); prof_pw_acc += uwTick - _t;
    /* Row 1 */
    _t = uwTick; poly_rej_ntt(poly1, rho, (uint8_t)(row0 + 1), (uint8_t)j); prof_rej_ntt += uwTick - _t;
    _t = uwTick; pointwise_acc_to(PKE_ACC1, poly1, poly0); prof_pw_acc += uwTick - _t;
    /* Row 2 */
    _t = uwTick; poly_rej_ntt(poly1, rho, (uint8_t)(row0 + 2), (uint8_t)j); prof_rej_ntt += uwTick - _t;
    _t = uwTick; pointwise_acc_to(PKE_ACC2, poly1, poly0); prof_pw_acc += uwTick - _t;
  }
}

/* ------------------------------------------------------------------ */
/*  ML-DSA-65 Signing (FIPS 204, Algorithm 7)                         */
/*                                                                    */
/*  Three-pass approach to minimize memory:                           */
/*    Pass 1: Compute w = Ay, hash w1 → c_tilde.                      */
/*    Pass 2: Compute z = y + c*s1, check ||z||∞ < γ1 − β.            */
/*    Pass 3: Recompute w, check r0 = LowBits(w − c*s2),              */
/*            compute c*t0, make hint.                                */
/* ------------------------------------------------------------------ */

int ml_dsa_65_sign(uint8_t *sig, size_t *sig_len,
                   const uint8_t *msg, size_t msg_len,
                   const uint8_t *ctx, size_t ctx_len,
                   const uint8_t *sk) {
  /* ---- Stack-allocated working memory ---- */
  poly poly0, poly1;                     /* 2 × 1024 = 2048 bytes */
  SHA3_CTX_T shake_ctx;                 /* ~208 bytes */
  uint8_t mu[64];
  uint8_t rho_prime[64];
  uint8_t c_tilde[MLDSA_C_TILDE_BYTES]; /* 48 bytes */
  uint16_t kappa;
  int reject;
  /* PKE register layout (no Kronecker in sign path → all 48 regs free):
   *   Pass 1: PKE_ACC0/1/2 — batched w_hat; last batch caches w[3..5]
   *   Pass 2: PKE-free (reordered to preserve cached w[3..5])
   *   Pass 3: Phase A uses cached w[3..5]; Phase B recomputes w[0..2]
   * Challenge c stored as sparse (pos+sign, 98 B) on stack.
   * Total stack ≈ 2048 + 208 + 64 + 64 + 48 + 98 + 55 + misc ≈ 2.7 KB */

  if (ctx_len > 255) return -1;

  /* ---- Step 1: Compute mu = H(tr || 0x00 || ctx_len || ctx || msg) ---- */
  {
    uint8_t hdr[2];
    shake256_init(&shake_ctx);
    shake_update(&shake_ctx, sk_tr(sk), MLDSA_TRBYTES);
    hdr[0] = 0x00;
    hdr[1] = (uint8_t)ctx_len;
    shake_update(&shake_ctx, hdr, 2);
    if (ctx_len > 0)
      shake_update(&shake_ctx, ctx, ctx_len);
    shake_update(&shake_ctx, msg, msg_len);
    shake_finalize(&shake_ctx);
    shake_squeeze(&shake_ctx, mu, 64);
  }

  /* ---- Step 2: Compute rho' = H(K || rnd || mu) ---- */
  {
    uint8_t rnd[32];
    memset(rnd, 0, 32);   /* deterministic signing */
    shake256_init(&shake_ctx);
    shake_update(&shake_ctx, sk_K(sk), 32);
    shake_update(&shake_ctx, rnd, 32);
    shake_update(&shake_ctx, mu, 64);
    shake_finalize(&shake_ctx);
    shake_squeeze(&shake_ctx, rho_prime, 64);
  }

  /* ---- Step 3: Signing loop ---- */
  kappa = 0;
  int iteration = 0;
  for (;;) {
    reject = 0;
    uint32_t t_iter_start = uwTick;
    prof_expand_mask = prof_ntt_fwd = prof_rej_ntt = prof_pw_acc = prof_invntt = 0;

    /* ==== Pass 1: Compute c_tilde = H(mu || w1Encode(w)) ==== */
    uint32_t t_pass1_start = uwTick;
    shake256_init(&shake_ctx);
    shake_update(&shake_ctx, mu, 64);

    {
      static const int accs[3] = {PKE_ACC0, PKE_ACC1, PKE_ACC2};
      for (int batch = 0; batch < K; batch += 3) {
        /* Compute w_hat for 3 rows at once, y[j] expanded only once */
        compute_w_hat_batch3(sk_rho(sk), rho_prime, kappa, batch,
                             poly0, poly1);
        /* Extract each row, INTT, encode w1 */
        for (int b = 0; b < 3; b++) {
          pke_load_poly(poly0, accs[b]);
          { uint32_t _t = uwTick; invntt(poly0); prof_invntt += uwTick - _t; }
          poly_caddq(poly0);
          w1_encode_update(&shake_ctx, poly0);
          /* Last batch: cache w[i] back to PKE for Pass 3 Phase A */
          if (batch == K - 3)
            pke_store_poly(accs[b], poly0);
        }
      }
    }

    shake_finalize(&shake_ctx);
    shake_squeeze(&shake_ctx, c_tilde, MLDSA_C_TILDE_BYTES);

    /* ==== Compute challenge c (sparse) ==== */
    challenge_sparse_t ch;
    poly_challenge(poly0, c_tilde);
    challenge_to_sparse(&ch, poly0);
    uint32_t t_pass1_end = uwTick;

    /* ==== Pass 2: Compute z = y + c*s1, check bounds ==== */
    uint32_t t_pass2_start = uwTick;
    uint32_t t_pass2_mul = 0;
    /* Also encode z into signature as we go.
     * Reordered: compute c*s1_j first (into poly1), then expand y_j,
     * so PKE registers (with cached w[3..5]) are never touched. */
    for (int j = 0; j < L; j++) {
      /* s1_j → poly0 (centered [-ETA,ETA]), c*s1_j → poly1 (centered [-196,196]) */
      unpack_eta(poly0, sk_s1(sk, j));
      { uint32_t _tm0 = uwTick;
      poly_mul_sparse_centered(poly1, poly0, &ch);
      t_pass2_mul += uwTick - _tm0; }
      /* y_j → poly0 (centered [-(GAMMA1-1), GAMMA1]) */
      poly_expand_mask(poly0, rho_prime, kappa + (uint16_t)j);
      /* z_j = y_j + c*s1_j (both centered — no reduction needed) */
      poly_add(poly0, poly0, poly1);
      /* Check ||z_j||∞ >= gamma1 - beta → reject */
      for (int n = 0; n < N; n++) {
        if (poly0[n] >= GAMMA1 - BETA_B || poly0[n] <= -(GAMMA1 - BETA_B)) {
          reject = 1;
          break;
        }
      }
      if (reject) break;
      /* Encode z_j (now centered) directly into signature */
      pack_z(sig + MLDSA_C_TILDE_BYTES + j * MLDSA_POLYZ_PACKEDBYTES,
             poly0);
    }
    if (reject) {
      uint32_t t_pass2_end = uwTick;
      printf("[PROF] iter=%d REJECTED@pass2 pass1=%lu pass2=%lu (mul=%lu) [em=%lu ntt=%lu rej=%lu pw=%lu intt=%lu]\n",
             iteration, (unsigned long)(t_pass1_end - t_pass1_start),
             (unsigned long)(t_pass2_end - t_pass2_start),
             (unsigned long)t_pass2_mul,
             (unsigned long)prof_expand_mask, (unsigned long)prof_ntt_fwd,
             (unsigned long)prof_rej_ntt, (unsigned long)prof_pw_acc,
             (unsigned long)prof_invntt);
      iteration++; kappa += L; continue;
    }
    uint32_t t_pass2_end = uwTick;

    /* ==== Pass 3: Check r0 & ct0, make hint ==== */
    /* Rows 3-5: use cached w from PKE (saved in Pass 1 last batch).
     * Rows 0-2: recompute w via compute_w_hat_batch3.
     * Process 3-5 first (while cache valid), then 0-2.
     * Hints stored out-of-order and merged at the end. */
    uint32_t t_pass3_start = uwTick;
    uint32_t t_pass3_w = 0;    /* time in compute_w_hat_row + INTT */
    uint32_t t_pass3_mul = 0;  /* time in poly_mul_sparse */
    {
      uint8_t *hint_buf = sig + MLDSA_C_TILDE_BYTES
                        + L * MLDSA_POLYZ_PACKEDBYTES;
      memset(hint_buf, 0, OMEGA + K);

      /* Temp hint storage for rows 3-5 (processed first from cache) */
      uint8_t temp_hints[OMEGA];
      int temp_count = 0;
      int temp_row_end[3];

      static const int accs[3] = {PKE_ACC0, PKE_ACC1, PKE_ACC2};

      /* --- Phase A: rows 3-5, w cached in PKE_ACC0/1/2 from Pass 1 --- */
      for (int b = 0; b < 3 && !reject; b++) {
        int i = 3 + b;
        /* w[i] already in accs[b] from Pass 1 cache — no recompute.
         * Compute c*s2[i]: s2 → poly0 (centered), result → poly1 (centered) */
        unpack_eta(poly0, sk_s2(sk, i));
        { uint32_t _tm0 = uwTick;
        poly_mul_sparse_centered(poly1, poly0, &ch);
        t_pass3_mul += uwTick - _tm0; }

        /* r = w[i] - c*s2[i] */
        pke_load_poly(poly0, accs[b]);
        poly_sub(poly0, poly0, poly1);
        poly_reduce(poly0);

        /* Compute c*t0[i] before the check loop so we can fuse
         * LowBits check + ct0 bounds + hint into a single pass,
         * doing decompose(r) only once per coefficient.
         * poly_mul_sparse_centered is < 1 ms — negligible cost. */
        pke_store_poly(accs[b], poly0);
        unpack_t0(poly0, sk_t0(sk, i));
        { uint32_t _tm0 = uwTick;
        poly_mul_sparse_centered(poly1, poly0, &ch);
        t_pass3_mul += uwTick - _tm0; }
        pke_load_poly(poly0, accs[b]);

        /* Fused loop: one decompose(r) per coeff gives r0 and r1.
         * r0 → LowBits reject, ct0 → bounds reject,
         * make_hint_from_decomposed(r1, r0, ct0) → pure-threshold hint
         * with zero additional decompose/freeze/high_bits calls. */
        for (int n = 0; n < N; n++) {
          int32_t r0;
          int32_t r1 = decompose(&r0, poly0[n]);
          if (r0 >= (int32_t)(GAMMA2 - BETA_B) ||
              r0 <= -(int32_t)(GAMMA2 - BETA_B)) {
            reject = 1;
            break;
          }
          if (poly1[n] >= (int32_t)GAMMA2 || poly1[n] <= -(int32_t)GAMMA2) {
            reject = 1;
            break;
          }
          if (make_hint_from_decomposed(r1, r0, poly1[n])) {
            if (temp_count >= OMEGA) { reject = 1; break; }
            temp_hints[temp_count++] = (uint8_t)n;
          }
        }
        if (reject) break;
        temp_row_end[b] = temp_count;
      }

      /* --- Phase B: rows 0-2, recompute w --- */
      int hint_count = 0;
      if (!reject) {
        { uint32_t _tw0 = uwTick;
        compute_w_hat_batch3(sk_rho(sk), rho_prime, kappa, 0,
                             poly0, poly1);
        t_pass3_w += uwTick - _tw0; }

        for (int b = 0; b < 3 && !reject; b++) {
          int i = b;
          /* Extract w_hat[i], INTT → w[i] */
          { uint32_t _tw0 = uwTick;
          pke_load_poly(poly0, accs[b]);
          { uint32_t _t = uwTick; invntt(poly0); prof_invntt += uwTick - _t; }
          poly_caddq(poly0);
          t_pass3_w += uwTick - _tw0; }
          pke_store_poly(accs[b], poly0);

          /* Compute c*s2[i]: s2 → poly0 (centered), result → poly1 (centered) */
          unpack_eta(poly0, sk_s2(sk, i));
          { uint32_t _tm0 = uwTick;
          poly_mul_sparse_centered(poly1, poly0, &ch);
          t_pass3_mul += uwTick - _tm0; }

          /* Restore w[i], compute r = w[i] - c*s2[i] */
          pke_load_poly(poly0, accs[b]);
          poly_sub(poly0, poly0, poly1);
          poly_reduce(poly0);

          /* Compute c*t0[i] before the check loop (same fuse as Phase A). */
          pke_store_poly(accs[b], poly0);
          unpack_t0(poly0, sk_t0(sk, i));
          { uint32_t _tm0 = uwTick;
          poly_mul_sparse_centered(poly1, poly0, &ch);
          t_pass3_mul += uwTick - _tm0; }
          pke_load_poly(poly0, accs[b]);

          /* Fused loop: one decompose(r) per coeff (same as Phase A). */
          for (int n = 0; n < N; n++) {
            int32_t r0;
            int32_t r1 = decompose(&r0, poly0[n]);
            if (r0 >= (int32_t)(GAMMA2 - BETA_B) ||
                r0 <= -(int32_t)(GAMMA2 - BETA_B)) {
              reject = 1;
              break;
            }
            if (poly1[n] >= (int32_t)GAMMA2 || poly1[n] <= -(int32_t)GAMMA2) {
              reject = 1;
              break;
            }
            if (make_hint_from_decomposed(r1, r0, poly1[n])) {
              if (hint_count + temp_count >= OMEGA) { reject = 1; break; }
              hint_buf[hint_count++] = (uint8_t)n;
            }
          }
          if (reject) break;
          hint_buf[OMEGA + i] = (uint8_t)hint_count;
        }
      }

      /* --- Merge: append rows 3-5 hints after rows 0-2 --- */
      if (!reject) {
        memcpy(hint_buf + hint_count, temp_hints, temp_count);
        hint_buf[OMEGA + 3] = (uint8_t)(hint_count + temp_row_end[0]);
        hint_buf[OMEGA + 4] = (uint8_t)(hint_count + temp_row_end[1]);
        hint_buf[OMEGA + 5] = (uint8_t)(hint_count + temp_row_end[2]);
      }
    }
    if (reject) {
      uint32_t t_pass3_end = uwTick;
      printf("[PROF] iter=%d REJECTED@pass3 pass1=%lu pass2=%lu (mul=%lu) pass3=%lu (w=%lu mul=%lu) [em=%lu ntt=%lu rej=%lu pw=%lu intt=%lu]\n",
             iteration,
             (unsigned long)(t_pass1_end - t_pass1_start),
             (unsigned long)(t_pass2_end - t_pass2_start),
             (unsigned long)t_pass2_mul,
             (unsigned long)(t_pass3_end - t_pass3_start),
             (unsigned long)t_pass3_w, (unsigned long)t_pass3_mul,
             (unsigned long)prof_expand_mask, (unsigned long)prof_ntt_fwd,
             (unsigned long)prof_rej_ntt, (unsigned long)prof_pw_acc,
             (unsigned long)prof_invntt);
      iteration++; kappa += L; continue;
    }
    {
      uint32_t t_pass3_end = uwTick;
      printf("[PROF] iter=%d ACCEPTED pass1=%lu pass2=%lu (mul=%lu) pass3=%lu (w=%lu mul=%lu) total=%lu [em=%lu ntt=%lu rej=%lu pw=%lu intt=%lu]\n",
             iteration,
             (unsigned long)(t_pass1_end - t_pass1_start),
             (unsigned long)(t_pass2_end - t_pass2_start),
             (unsigned long)t_pass2_mul,
             (unsigned long)(t_pass3_end - t_pass3_start),
             (unsigned long)t_pass3_w, (unsigned long)t_pass3_mul,
             (unsigned long)(t_pass3_end - t_iter_start),
             (unsigned long)prof_expand_mask, (unsigned long)prof_ntt_fwd,
             (unsigned long)prof_rej_ntt, (unsigned long)prof_pw_acc,
             (unsigned long)prof_invntt);
    }

    /* ==== Encode signature ==== */
    /* c_tilde is at the beginning */
    memcpy(sig, c_tilde, MLDSA_C_TILDE_BYTES);
    /* z is already encoded in pass 2 */
    /* hint is already encoded in pass 3 */

    *sig_len = MLDSA_SIG_BYTES;
    return 0;
  }
}

/* ------------------------------------------------------------------ */
/*  ML-DSA-65 Signing from seed (no sk buffer)                        */
/*                                                                    */
/*  Same three-pass approach as ml_dsa_65_sign, but regenerates s1,   */
/*  s2, t0 from seed on the fly.  ~2-3x slower due to recomputation.  */
/* ------------------------------------------------------------------ */

/* Helper: regenerate t0[i] from keygen secrets.
 * Computes row i of t = A*s1 + s2, then Power2Round to extract t0.
 * WARNING: clobbers all PKE registers.
 * Result in out[N], coefficients in [0, Q).
 * scratch[N] is used as temporary space. */
static void regen_t0(int32_t out[N], const uint8_t rho[32],
                     const uint8_t rho_prime[64], int i) {
  /* Accumulate A_hat[i][j] * s1_hat[j] in PKE_SLOT1.
   * Uses only out[N] as scratch — s1_hat stored in PKE_SLOT0,
   * read by pointwise_acc_pke_slot0. */
  {
    uint32_t zeros[PKE_COEFFS_PER_REG];
    memset(zeros, 0, sizeof(zeros));
    for (int r = 0; r < N / PKE_COEFFS_PER_REG; r++)
      eccPKEWriteBuf(PKE_SLOT1 + r, zeros, PKE_COEFFS_PER_REG);
  }
  for (int j = 0; j < L; j++) {
    poly_rej_bounded(out, rho_prime, (uint16_t)j);
    ntt(out);
    pke_store_poly(PKE_SLOT0, out);            /* s1_hat → SLOT0 */
    poly_rej_ntt(out, rho, (uint8_t)i, (uint8_t)j);  /* A_hat → out */
    pointwise_acc_pke_slot0(out);              /* SLOT1 += A_hat * SLOT0 */
  }
  pke_load_poly(out, PKE_SLOT1);
  invntt(out);

  /* Add s2[i]: save t to SLOT1, gen s2 into out, add from SLOT1 */
  pke_store_poly(PKE_SLOT1, out);
  poly_rej_bounded(out, rho_prime, (uint16_t)(L + i));
  {
    int32_t t_chunk[PKE_COEFFS_PER_REG];
    for (int r = 0; r < N / PKE_COEFFS_PER_REG; r++) {
      eccPKEReadBuf((uint32_t *)t_chunk, PKE_SLOT1 + r, PKE_COEFFS_PER_REG);
      for (int k = 0; k < PKE_COEFFS_PER_REG; k++)
        out[r * PKE_COEFFS_PER_REG + k] += t_chunk[k];
    }
  }
  poly_caddq(out);

  /* Power2Round: extract t0 */
  for (int n = 0; n < N; n++) {
    int32_t t0_coeff;
    (void)power2round(&t0_coeff, out[n]);
    out[n] = t0_coeff;
  }
  poly_caddq(out);
}

/* Helper: restore challenge c from c_tilde into PKE_SLOT0.
 * Uses scratch[N] as temporary. */
static void restore_challenge(const uint8_t c_tilde[MLDSA_C_TILDE_BYTES],
                              int32_t scratch[N]) {
  poly_challenge(scratch, c_tilde);
  poly_caddq(scratch);
  pke_store_poly(PKE_SLOT0, scratch);
}

/* ------------------------------------------------------------------ */

/* Re-derive rho and rho_prime_keygen from seed.
 * Marked noinline so the ~200B SHA3 context is freed on return,
 * rather than inflating the caller's stack frame. */
__attribute__((noinline))
static void derive_keygen_secrets(const uint8_t seed[32],
                                  uint8_t rho[32],
                                  uint8_t rho_prime_keygen[64],
                                  uint8_t K_seed[32]) {
  SHA3_CTX_T ctx;
  uint8_t dom[2] = { (uint8_t)K, (uint8_t)L };
  shake256_init(&ctx);
  shake_update(&ctx, seed, 32);
  shake_update(&ctx, dom, 2);
  shake_finalize(&ctx);
  shake_squeeze(&ctx, rho, 32);
  shake_squeeze(&ctx, rho_prime_keygen, 64);
  if (K_seed) shake_squeeze(&ctx, K_seed, 32);
}

int ml_dsa_65_sign_seed(uint8_t *sig, size_t *sig_len,
                        const uint8_t *msg, size_t msg_len,
                        const uint8_t *ctx, size_t ctx_len,
                        const uint8_t *seed, const uint8_t *tr) {
  poly poly0, poly1;
  SHA3_CTX_T shake_ctx;
  uint8_t mu[64];
  uint8_t rho[32], rho_prime_sign[64], K_seed[32];
  uint8_t rho_prime_keygen[64];
  uint8_t c_tilde[MLDSA_C_TILDE_BYTES];
  uint16_t kappa;
  int reject;

  if (ctx_len > 255) return -1;

  /* ---- Derive keygen secrets from seed ---- */
  derive_keygen_secrets(seed, rho, rho_prime_keygen, K_seed);

  /* ---- Compute mu = H(tr || 0x00 || ctx_len || ctx || msg) ---- */
  {
    uint8_t hdr[2];
    shake256_init(&shake_ctx);
    shake_update(&shake_ctx, tr, MLDSA_TRBYTES);
    hdr[0] = 0x00;
    hdr[1] = (uint8_t)ctx_len;
    shake_update(&shake_ctx, hdr, 2);
    if (ctx_len > 0)
      shake_update(&shake_ctx, ctx, ctx_len);
    shake_update(&shake_ctx, msg, msg_len);
    shake_finalize(&shake_ctx);
    shake_squeeze(&shake_ctx, mu, 64);
  }

  /* ---- Compute rho' for signing = H(K || rnd || mu) ---- */
  {
    uint8_t rnd[32];
    memset(rnd, 0, 32);
    shake256_init(&shake_ctx);
    shake_update(&shake_ctx, K_seed, 32);
    shake_update(&shake_ctx, rnd, 32);
    shake_update(&shake_ctx, mu, 64);
    shake_finalize(&shake_ctx);
    shake_squeeze(&shake_ctx, rho_prime_sign, 64);
  }

  /* ---- Signing loop ---- */
  kappa = 0;
  for (;;) {
    reject = 0;

    /* ==== Pass 1: c_tilde = H(mu || w1Encode(w)) ==== */
    shake256_init(&shake_ctx);
    shake_update(&shake_ctx, mu, 64);
    for (int i = 0; i < K; i++) {
      compute_w_hat_row(rho, rho_prime_sign, kappa, i, poly0, poly1);
      pke_load_poly(poly0, PKE_SLOT1);
      invntt(poly0);
      poly_caddq(poly0);
      w1_encode_update(&shake_ctx, poly0);
    }
    shake_finalize(&shake_ctx);
    shake_squeeze(&shake_ctx, c_tilde, MLDSA_C_TILDE_BYTES);

    /* ==== Challenge c ==== */
    restore_challenge(c_tilde, poly0);

    /* ==== Pass 2: z = y + c*s1, check bounds ==== */
    for (int j = 0; j < L; j++) {
      poly_expand_mask(poly0, rho_prime_sign, kappa + (uint16_t)j);
      pke_store_poly(PKE_SLOT1, poly0);
      poly_rej_bounded(poly0, rho_prime_keygen, (uint16_t)j);
      poly_caddq(poly0);
      poly_mul_challenge(poly1, poly0);
      pke_load_poly(poly0, PKE_SLOT1);
      poly_add(poly0, poly0, poly1);
      poly_reduce(poly0);
      for (int n = 0; n < N; n++) {
        if (poly0[n] > Q / 2) poly0[n] -= Q;
        if (poly0[n] >= GAMMA1 - BETA_B || poly0[n] <= -(GAMMA1 - BETA_B)) {
          reject = 1;
          break;
        }
      }
      if (reject) break;
      pack_z(sig + MLDSA_C_TILDE_BYTES + j * MLDSA_POLYZ_PACKEDBYTES,
             poly0);
    }
    if (reject) { kappa += L; continue; }

    /* ==== Pass 3: check r0, ct0, make hint ==== */
    /* Per row:
     *   Phase A: compute r = w-cs2, check LowBits (uses poly0+poly1)
     *   Phase B: regen t0 (single-scratch), compute ct0 → poly1,
     *     check ct0 bounds, recompute w using pointwise_acc_pke_slot0
     *     (preserves poly1=ct0), compute cs2 into stack tmp, get r,
     *     make hint. */
    {
      uint8_t *hint_buf = sig + MLDSA_C_TILDE_BYTES
                        + L * MLDSA_POLYZ_PACKEDBYTES;
      int hint_count = 0;
      memset(hint_buf, 0, OMEGA + K);

      for (int i = 0; i < K; i++) {
        /* Phase A: r = w - cs2, LowBits check */
        compute_w_hat_row(rho, rho_prime_sign, kappa, i, poly0, poly1);
        pke_load_poly(poly0, PKE_SLOT1);
        invntt(poly0);
        poly_caddq(poly0);
        pke_store_poly(PKE_SLOT1, poly0);  /* save w */

        restore_challenge(c_tilde, poly0);
        poly_rej_bounded(poly0, rho_prime_keygen, (uint16_t)(L + i));
        poly_caddq(poly0);
        poly_mul_challenge(poly1, poly0);  /* cs2 */

        pke_load_poly(poly0, PKE_SLOT1);   /* w */
        poly_sub(poly0, poly0, poly1);
        poly_reduce(poly0);

        for (int n = 0; n < N; n++) {
          int32_t r0;
          decompose(&r0, poly0[n]);
          if (r0 >= (int32_t)(GAMMA2 - BETA_B) ||
              r0 <= -(int32_t)(GAMMA2 - BETA_B)) {
            reject = 1;
            break;
          }
        }
        if (reject) break;

        /* Phase B: ct0 + hint */
        regen_t0(poly0, rho, rho_prime_keygen, i);

        restore_challenge(c_tilde, poly1);
        poly_mul_challenge(poly1, poly0);  /* ct0 — keep in poly1 */

        for (int n = 0; n < N; n++) {
          int32_t v = poly1[n];
          if (v > Q / 2) v -= Q;
          if (v >= (int32_t)GAMMA2 || v <= -(int32_t)GAMMA2) {
            reject = 1;
            break;
          }
        }
        if (reject) break;

        /* Recompute w using pointwise_acc_pke_slot0 (only poly0,
         * preserves poly1=ct0) */
        {
          uint32_t zeros[PKE_COEFFS_PER_REG];
          memset(zeros, 0, sizeof(zeros));
          for (int r = 0; r < N / PKE_COEFFS_PER_REG; r++)
            eccPKEWriteBuf(PKE_SLOT1 + r, zeros, PKE_COEFFS_PER_REG);
        }
        for (int j = 0; j < L; j++) {
          poly_expand_mask(poly0, rho_prime_sign, kappa + (uint16_t)j);
          ntt(poly0);
          pke_store_poly(PKE_SLOT0, poly0);
          poly_rej_ntt(poly0, rho, (uint8_t)i, (uint8_t)j);
          pointwise_acc_pke_slot0(poly0);
        }
        pke_load_poly(poly0, PKE_SLOT1);
        invntt(poly0);
        poly_caddq(poly0);

        /* cs2 → stack tmp (avoids aliased poly_mul_challenge) */
        pke_store_poly(PKE_SLOT1, poly0);  /* save w */
        restore_challenge(c_tilde, poly0);
        poly_rej_bounded(poly0, rho_prime_keygen, (uint16_t)(L + i));
        poly_caddq(poly0);
        {
          int32_t cs2[N];
          poly_mul_challenge(cs2, poly0);
          pke_load_poly(poly0, PKE_SLOT1);  /* w */
          poly_sub(poly0, poly0, cs2);
        }
        poly_reduce(poly0);
        /* poly0 = r, poly1 = ct0 */

        /* One decompose(r) per coeff; pure-threshold hint.
         * poly1 = ct0 in [0,Q) from poly_mul_challenge — center it. */
        for (int n = 0; n < N; n++) {
          int32_t r0;
          int32_t r1 = decompose(&r0, poly0[n]);
          int32_t ct0 = poly1[n];
          if (ct0 > Q / 2) ct0 -= Q;
          if (make_hint_from_decomposed(r1, r0, ct0)) {
            if (hint_count >= OMEGA) { reject = 1; break; }
            hint_buf[hint_count++] = (uint8_t)n;
          }
        }
        if (reject) break;
        hint_buf[OMEGA + i] = (uint8_t)hint_count;
      }
    }
    if (reject) { kappa += L; continue; }

    /* ==== Encode signature ==== */
    memcpy(sig, c_tilde, MLDSA_C_TILDE_BYTES);
    *sig_len = MLDSA_SIG_BYTES;
    return 0;
  }
}

/*                                                                    */
/*  pk = (rho || t1_packed)                                           */
/*  sk = (rho || K || tr || s1_packed || s2_packed || t0_packed)      */
/*                                                                    */
/*  Memory: reuses 3 stack polys + PKE buffer.                        */
/*  t = NTT^{-1}(A_hat * NTT(s1)) + s2, then Power2Round.             */
/* ------------------------------------------------------------------ */

int ml_dsa_65_keygen(uint8_t *pk, uint8_t *sk, uint8_t *tr_out,
                     const uint8_t *seed) {
  poly poly0, poly1, poly2;
  SHA3_CTX_T tr_ctx;       /* streaming H(pk) for tr */
  uint8_t rho[32], rho_prime[64], K_seed[32];
  int compute_tr = (tr_out != NULL || (sk != NULL));  /* sk needs tr too */

  /* Step 1: (rho, rho', K) = H(seed) */
  derive_keygen_secrets(seed, rho, rho_prime, K_seed);

  /* Copy rho and K into sk */
  if (sk) {
    memcpy(sk, rho, 32);           /* sk[0..31]  = rho */
    memcpy(sk + 32, K_seed, 32);   /* sk[32..63] = K   */
  }

  /* Step 2: Generate s1, s2 using rho' and pack into sk */
  for (int j = 0; j < L; j++) {
    poly_rej_bounded(poly0, rho_prime, (uint16_t)j);
    if (sk)
      pack_eta(sk + 128 + j * MLDSA_POLYETA_PACKEDBYTES, poly0);
  }
  for (int i = 0; i < K; i++) {
    poly_rej_bounded(poly0, rho_prime, (uint16_t)(L + i));
    if (sk)
      pack_eta(sk + 128 + L * MLDSA_POLYETA_PACKEDBYTES
               + i * MLDSA_POLYETA_PACKEDBYTES, poly0);
  }

  /* Copy rho to pk, start tr hash */
  if (pk) memcpy(pk, rho, 32);
  if (compute_tr) {
    shake256_init(&tr_ctx);
    shake_update(&tr_ctx, rho, 32);
  }

  /* Step 3: Compute t = A*s1 + s2, row by row.
   * For each row i, compute t[i], then Power2Round → (t1, t0).
   * Pack t1 into pk (if non-NULL) and feed into tr hash.
   * Pack t0 into sk (if non-NULL). */
  for (int i = 0; i < K; i++) {
    /* t_hat[i] = sum_j A_hat[i][j] * s1_hat[j] */
    memset(poly2, 0, N * sizeof(int32_t));
    for (int j = 0; j < L; j++) {
      /* Regenerate s1[j] from rho' (sk may be NULL) */
      poly_rej_bounded(poly0, rho_prime, (uint16_t)j);
      ntt(poly0);
      poly_rej_ntt(poly1, rho, (uint8_t)i, (uint8_t)j);
      poly_pointwise_acc(poly2, poly1, poly0);
    }
    invntt(poly2);

    /* Add s2[i] */
    poly_rej_bounded(poly0, rho_prime, (uint16_t)(L + i));
    poly_add(poly2, poly2, poly0);
    poly_caddq(poly2);

    /* Power2Round: t = t1 * 2^D + t0 */
    for (int n = 0; n < N; n++) {
      int32_t t0_coeff;
      poly0[n] = power2round(&t0_coeff, poly2[n]); /* poly0 = t1 */
      poly1[n] = t0_coeff;                          /* poly1 = t0 */
    }

    /* Pack t1 into pk and/or feed into tr hash */
    {
      uint8_t t1_packed[MLDSA_POLYT1_PACKEDBYTES];
      pack_t1(t1_packed, poly0);
      if (pk)
        memcpy(pk + 32 + i * MLDSA_POLYT1_PACKEDBYTES,
               t1_packed, MLDSA_POLYT1_PACKEDBYTES);
      if (compute_tr)
        shake_update(&tr_ctx, t1_packed, MLDSA_POLYT1_PACKEDBYTES);
    }

    /* Pack t0 into sk */
    if (sk)
      pack_t0(sk + 128 + (L + K) * MLDSA_POLYETA_PACKEDBYTES
              + i * MLDSA_POLYT0_PACKEDBYTES, poly1);
  }

  /* Step 4: Finalize tr = H(pk) */
  if (compute_tr) {
    uint8_t tr_buf[MLDSA_TRBYTES];
    shake_finalize(&tr_ctx);
    shake_squeeze(&tr_ctx, tr_buf, MLDSA_TRBYTES);
    if (tr_out)
      memcpy(tr_out, tr_buf, MLDSA_TRBYTES);
    if (sk)
      memcpy(sk + 64, tr_buf, MLDSA_TRBYTES);
  }

  return 0;
}

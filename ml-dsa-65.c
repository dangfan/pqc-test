/*
 * ML-DSA-65 (FIPS 204) - Sign-only implementation.
 * Memory budget: <= 6 KB. Matrix A expanded on-the-fly.
 * Uses hardware MMM (sm2MonMul) and PKE registers as data buffer.
 */

#include "ml-dsa-65.h"
#include "se.h"
#include "sha3.h"
#include <string.h>
#include <stdio.h>

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
/*  Modular arithmetic helpers                                         */
/* ------------------------------------------------------------------ */

static int32_t freeze(int32_t a) {
  a %= Q;
  if (a < 0) a += Q;
  return a;
}

/* ------------------------------------------------------------------ */
/*  PKE register read/write for polynomial storage                     */
/*  Store coefficients as native uint32_t words (16 per register).     */
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

/* ------------------------------------------------------------------ */
/*  Host ↔ PKE byte-order helpers                                     */
/*                                                                      */
/*  PKE registers store multi-precision integers in big-endian byte     */
/*  order: the most significant byte is at the lowest address.          */
/*  On the host side, scalar uint32_t values use native byte order,     */
/*  and multi-word integers use word[0] = least significant (LE word    */
/*  order).  The helpers below convert between the two conventions.     */
/* ------------------------------------------------------------------ */

/* Convert a single uint32_t between host byte order and big-endian. */
static inline uint32_t host_to_be32(uint32_t x) {
#if __BYTE_ORDER__ == __ORDER_BIG_ENDIAN__
  return x;
#else
  return __builtin_bswap32(x);
#endif
}
#define be32_to_host(x) host_to_be32(x)  /* same operation, symmetric */

/* Write 32-bit value into upper half of 8-byte PKE register (value * 2^32). */
static void write_pke_hi32(int reg_idx, uint32_t val) {
  uint32_t buf[2] = { host_to_be32(val), 0 };
  eccPKEWriteBuf(reg_idx, buf, 2);
}

/* Write 32-bit value into lower half of 8-byte PKE register. */
static void write_pke_lo32(int reg_idx, uint32_t val) {
  uint32_t buf[2] = { 0, host_to_be32(val) };
  eccPKEWriteBuf(reg_idx, buf, 2);
}

/* Read 32-bit value from lower half of 8-byte PKE register. */
static uint32_t read_pke_lo32(int reg_idx) {
  uint32_t buf[2];
  eccPKEReadBuf(buf, reg_idx, 2);
  return be32_to_host(buf[1]);
}

static void pke_mmm_setup(void) {
  rsaPKESetLen(2);
  write_pke_lo32(0, Q);
  uint32_t mod_words[2] = { 0, host_to_be32(Q) };
  uint32_t mc[2];
  eccCalMc64(mc, mod_words);
  rsaPKEWriteMc(mc);
}

/* Hardware modular multiply: a * b * R_32^{-1} mod Q, where R_32 = 2^32.
 * With pkeLen=2, hardware R = 2^64. We store operand a in the UPPER half
 * (value = a * 2^32) and b in the LOWER half (value = b), so:
 *   MonMul(a*2^32, b) = a*2^32 * b * 2^{-64} mod Q = a * b * 2^{-32} mod Q
 * This matches the original R=2^32 Montgomery behavior. */
static int32_t hw_mon_mul(int32_t a, int32_t b) {
  write_pke_hi32(1, (uint32_t)freeze(a));
  write_pke_lo32(2, (uint32_t)freeze(b));
  sm2MonMul(3, 1, 2);
  return (int32_t)read_pke_lo32(3);
}

/* ------------------------------------------------------------------ */
/*  NTT / INTT (Cooley-Tukey / Gentleman-Sande, using hw MMM)         */
/* ------------------------------------------------------------------ */

static void ntt(int32_t a[N]) {
  unsigned int len, start, j, k;
  int32_t t;

  pke_mmm_setup();
  k = 0;
  for (len = 128; len >= 1; len >>= 1) {
    for (start = 0; start < N; start = j + len) {
      int32_t zeta = zetas[++k];
      for (j = start; j < start + len; ++j) {
        t = hw_mon_mul(zeta, a[j + len]);
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
    a[j] %= Q;

  pke_mmm_setup();
  k = 256;
  for (len = 1; len <= 128; len <<= 1) {
    for (start = 0; start < N; start = j + len) {
      int32_t zeta = -zetas[--k];
      for (j = start; j < start + len; ++j) {
        t = a[j];
        a[j]       = t + a[j + len];
        a[j + len] = t - a[j + len];
        a[j + len] = hw_mon_mul(zeta, a[j + len]);
      }
    }
  }
  for (j = 0; j < N; ++j)
    a[j] = hw_mon_mul(INTT_F, a[j]);
}

/* Pointwise multiplication: c = a * b (NTT domain, Montgomery) */
static void poly_pointwise(int32_t c[N], const int32_t a[N],
                           const int32_t b[N]) {
  pke_mmm_setup();
  for (int i = 0; i < N; i++)
    c[i] = hw_mon_mul(a[i], b[i]);
}

/* c += a * b (NTT domain, Montgomery, accumulate) */
static void poly_pointwise_acc(int32_t c[N], const int32_t a[N],
                               const int32_t b[N]) {
  pke_mmm_setup();
  for (int i = 0; i < N; i++)
    c[i] += hw_mon_mul(a[i], b[i]);
}

/* ------------------------------------------------------------------ */
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
/*    - Total: 256 HW multiplies (vs ~2304 in NTT approach).           */
/* ------------------------------------------------------------------ */

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

/* Unpack (2*KRON_G - 1) coefficients from a product big integer.
 * Each coefficient is at KRON_L-bit stride. Coefficients fit in 50 bits. */
static void kron_unpack(int64_t *coeffs, int ncoeffs,
                        const uint32_t prod[KRON_PROD_WORDS]) {
  uint64_t mask = ((uint64_t)1 << KRON_L) - 1;
  for (int i = 0; i < ncoeffs; i++) {
    int bit_pos = i * KRON_L;
    int word = bit_pos / 32;
    int shift = bit_pos % 32;
    uint64_t val = (uint64_t)prod[word] >> shift;
    if (shift + KRON_L > 32 && word + 1 < KRON_PROD_WORDS)
      val |= (uint64_t)prod[word + 1] << (32 - shift);
    if (shift + KRON_L > 64 && word + 2 < KRON_PROD_WORDS)
      val |= (uint64_t)prod[word + 2] << (64 - shift);
    coeffs[i] = (int64_t)(val & mask);
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

static void kron_bigint_mul(uint32_t result[KRON_PROD_WORDS],
                            const uint32_t a[KRON_PACK_WORDS],
                            const uint32_t b[KRON_PACK_WORDS]) {
  uint32_t tmp[KRON_PKE_LEN];

  /* Write a to PKE reg (host → PKE, zero-padded) */
  memset(tmp, 0, KRON_PKE_LEN * 4);
  memcpy(tmp, a, KRON_PACK_WORDS * 4);
  host_to_pke(tmp, KRON_PKE_LEN);
  eccPKEWriteBuf(KRON_REG_A, tmp, KRON_PKE_LEN);

  /* Write b to PKE reg (host → PKE, zero-padded) */
  memset(tmp, 0, KRON_PKE_LEN * 4);
  memcpy(tmp, b, KRON_PACK_WORDS * 4);
  host_to_pke(tmp, KRON_PKE_LEN);
  eccPKEWriteBuf(KRON_REG_B, tmp, KRON_PKE_LEN);

  sm2MonMul(KRON_REG_R, KRON_REG_A, KRON_REG_B);

  /* Read result and convert PKE → host */
  eccPKEReadBuf(tmp, KRON_REG_R, KRON_PKE_LEN);
  pke_to_host(tmp, KRON_PKE_LEN);
  memcpy(result, tmp, KRON_PROD_WORDS * 4);
}

/* Multiply two polynomials in Z_q[X]/(X^N+1) using Kronecker substitution.
 * Both a and b must have coefficients in [0, Q). Result c has coefficients
 * reduced mod Q in [0, Q). */
static void poly_mul(int32_t c[N], const int32_t a[N], const int32_t b[N]) {
  int64_t acc[N];  /* accumulator with inline negacyclic reduction, ~2 KB */
  memset(acc, 0, sizeof(acc));

  kron_mmm_setup();

  /* Schoolbook on groups: for each pair (gi, gj), compute sub-product
   * via Kronecker and accumulate with negacyclic folding (X^N = -1). */
  for (int gi = 0; gi < KRON_T; gi++) {
    const int32_t *a_seg = &a[gi * KRON_G];
    uint32_t a_packed[KRON_PACK_WORDS];
    kron_pack(a_packed, a_seg, KRON_G);

    for (int gj = 0; gj < KRON_T; gj++) {
      const int32_t *b_seg = &b[gj * KRON_G];
      uint32_t b_packed[KRON_PACK_WORDS];
      kron_pack(b_packed, b_seg, KRON_G);

      uint32_t prod[KRON_PROD_WORDS];
      kron_bigint_mul(prod, a_packed, b_packed);

      int64_t sub_coeffs[2 * KRON_G - 1];
      kron_unpack(sub_coeffs, 2 * KRON_G - 1, prod);

      int base = (gi + gj) * KRON_G;
      for (int k = 0; k < 2 * KRON_G - 1; k++) {
        int idx = base + k;
        if (idx < N)
          acc[idx] += sub_coeffs[k];
        else
          acc[idx - N] -= sub_coeffs[k];  /* X^N = -1 */
      }
    }
  }

  /* Reduce mod Q */
  for (int i = 0; i < N; i++) {
    int64_t v = acc[i] % Q;
    if (v < 0) v += Q;
    c[i] = (int32_t)v;
  }
}

/* Multiply and accumulate: c += a * b mod (X^N+1, Q) */
static void poly_mul_acc(int32_t c[N], const int32_t a[N], const int32_t b[N]) {
  poly tmp;
  poly_mul(tmp, a, b);
  for (int i = 0; i < N; i++) {
    int32_t v = c[i] + tmp[i];
    v %= Q;
    if (v < 0) v += Q;
    c[i] = v;
  }
}

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

/* ------------------------------------------------------------------ */
/*  Basic polynomial operations                                        */
/* ------------------------------------------------------------------ */

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
  for (int i = 0; i < N; i++) { if (a[i] < 0) a[i] += Q; }
}

/* ------------------------------------------------------------------ */
/*  Sampling: RejNTTPoly (ExpandA), ExpandMask, RejBoundedPoly,        */
/*            SampleInBall                                              */
/* ------------------------------------------------------------------ */

/* Rejection-sample one NTT-domain polynomial from SHAKE128 (ExpandA).
 * FIPS 204, Algorithm 30 (RejNTTPoly). */
static void poly_rej_ntt(int32_t a[N], const uint8_t rho[32],
                         uint8_t row, uint8_t col) {
  SHA3_CTX_T ctx;
  uint8_t buf[3];
  int ctr = 0;

  shake128_init(&ctx);
  shake_update(&ctx, rho, 32);
  buf[0] = col;
  buf[1] = row;
  shake_update(&ctx, buf, 2);
  shake_finalize(&ctx);

  while (ctr < N) {
    shake_squeeze(&ctx, buf, 3);
    uint32_t t = ((uint32_t)buf[0]) | ((uint32_t)buf[1] << 8) |
                 ((uint32_t)(buf[2] & 0x7F) << 16);
    if (t < (uint32_t)Q)
      a[ctr++] = (int32_t)t;
  }
}

/* Expand masking vector component y_idx from rho'.
 * FIPS 204, Algorithm 36 (ExpandMask). gamma1 = 2^19, 20-bit packing. */
static void poly_expand_mask(int32_t a[N], const uint8_t rho_prime[64],
                             uint16_t idx) {
  SHA3_CTX_T ctx;
  uint8_t buf[5];
  uint8_t hdr[2];

  shake256_init(&ctx);
  shake_update(&ctx, rho_prime, 64);
  hdr[0] = (uint8_t)(idx & 0xFF);
  hdr[1] = (uint8_t)(idx >> 8);
  shake_update(&ctx, hdr, 2);
  shake_finalize(&ctx);

  /* 20-bit packed: 4 coefficients per 10 bytes → squeeze 5 bytes, get 2 */
  for (int i = 0; i < N / 2; i++) {
    shake_squeeze(&ctx, buf, 5);
    uint32_t t0 = ((uint32_t)buf[0]) | ((uint32_t)buf[1] << 8) |
                  ((uint32_t)(buf[2] & 0x0F) << 16);
    uint32_t t1 = ((uint32_t)buf[2] >> 4) | ((uint32_t)buf[3] << 4) |
                  ((uint32_t)buf[4] << 12);
    a[2 * i]     = (int32_t)(GAMMA1 - t0);
    a[2 * i + 1] = (int32_t)(GAMMA1 - t1);
  }
}

/* Rejection-sample a bounded polynomial with coefficients in [-eta, eta].
 * FIPS 204, Algorithm 31 (RejBoundedPoly). eta=4: sample nibbles. */
static void poly_rej_bounded(int32_t a[N], const uint8_t seed[64],
                             uint16_t nonce) {
  SHA3_CTX_T ctx;
  uint8_t buf[1];
  uint8_t hdr[2];
  int ctr = 0;

  shake256_init(&ctx);
  shake_update(&ctx, seed, 64);
  hdr[0] = (uint8_t)(nonce & 0xFF);
  hdr[1] = (uint8_t)(nonce >> 8);
  shake_update(&ctx, hdr, 2);
  shake_finalize(&ctx);

  while (ctr < N) {
    shake_squeeze(&ctx, buf, 1);
    uint8_t z0 = buf[0] & 0x0F;
    uint8_t z1 = buf[0] >> 4;
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
  uint8_t buf[8];
  uint64_t signs;

  memset(c, 0, N * sizeof(int32_t));

  shake256_init(&ctx);
  shake_update(&ctx, c_tilde, MLDSA_C_TILDE_BYTES);
  shake_finalize(&ctx);
  shake_squeeze(&ctx, buf, 8);
  signs = 0;
  for (int i = 0; i < 8; i++)
    signs |= (uint64_t)buf[i] << (8 * i);

  for (unsigned int i = N - TAU; i < N; i++) {
    uint8_t b;
    do {
      shake_squeeze(&ctx, &b, 1);
    } while ((unsigned int)b > i);
    unsigned int j = (unsigned int)b;
    c[i] = c[j];
    c[j] = (signs & 1) ? -1 : 1;
    signs >>= 1;
  }
}

/* ------------------------------------------------------------------ */
/*  Packing / Unpacking                                                */
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

/* Encode w1 polynomial (4 bits per coeff for gamma2=(q-1)/32). */
static void pack_w1(uint8_t *buf, const int32_t a[N]) {
  for (int i = 0; i < N / 2; i++)
    buf[i] = (uint8_t)((uint32_t)a[2*i] | ((uint32_t)a[2*i+1] << 4));
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

/* ------------------------------------------------------------------ */
/*  Decompose, HighBits, LowBits, MakeHint, Power2Round              */
/* ------------------------------------------------------------------ */

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

static int32_t low_bits(int32_t a) {
  int32_t r0;
  decompose(&r0, a);
  return r0;
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

/* Pack a w1 polynomial (already computed, 4-bit coefficients) and feed
 * directly into a SHAKE context.  Same streaming approach as
 * w1_encode_update but skips the HighBits step. */
static void __attribute__((noinline)) w1_pack_update(SHA3_CTX_T *ctx,
                                                     const int32_t w1[N]) {
  uint8_t chunk[16];
  for (int i = 0; i < N; i += 32) {
    for (int j = 0; j < 16; j++)
      chunk[j] = (uint8_t)((uint32_t)w1[i + 2*j]
                          | ((uint32_t)w1[i + 2*j + 1] << 4));
    shake_update(ctx, chunk, 16);
  }
}

static int make_hint(int32_t z, int32_t r) {
  return (high_bits(r) != high_bits(freeze(r + z)));
}

/* Power2Round: a = a1 * 2^D + a0 */
static int32_t power2round(int32_t *a0, int32_t a) {
  int32_t a1 = (a + (1 << (D_BITS - 1)) - 1) >> D_BITS;
  *a0 = a - (a1 << D_BITS);
  return a1;
}

/* ------------------------------------------------------------------ */
/*  Secret key layout accessors                                        */
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
/*  One row at a time. y is re-expanded for each row.                  */
/*  Uses poly0 for y_hat_j, poly1 for A_hat_ij.                       */
/*  w_hat is accumulated in PKE_SLOT1 (no stack poly needed).          */
/*  Caller retrieves result via pke_load_poly(dst, PKE_SLOT1).         */
/* ------------------------------------------------------------------ */

/* Pointwise accumulate into PKE_SLOT1: w_hat += a * b (NTT domain).
 * Reads/writes w_hat from PKE_SLOT1 in 16-coeff chunks.
 * Only touches regs 0-3 (MMM) and 20-35 (PKE_SLOT1). */
static void pointwise_acc_pke(const int32_t a[N], const int32_t b[N]) {
  pke_mmm_setup();
  int32_t w_chunk[PKE_COEFFS_PER_REG];
  for (int r = 0; r < N / PKE_COEFFS_PER_REG; r++) {
    eccPKEReadBuf((uint32_t *)w_chunk, PKE_SLOT1 + r, PKE_COEFFS_PER_REG);
    for (int k = 0; k < PKE_COEFFS_PER_REG; k++)
      w_chunk[k] += hw_mon_mul(a[r * PKE_COEFFS_PER_REG + k],
                                b[r * PKE_COEFFS_PER_REG + k]);
    eccPKEWriteBuf(PKE_SLOT1 + r, (const uint32_t *)w_chunk,
                   PKE_COEFFS_PER_REG);
  }
}

static void compute_w_hat_row(const uint8_t *rho,
                              const uint8_t *rho_prime,
                              uint16_t kappa, int row,
                              int32_t poly0[N], int32_t poly1[N]) {
  /* Zero PKE_SLOT1 (w_hat accumulator) */
  {
    uint32_t zeros[PKE_COEFFS_PER_REG];
    memset(zeros, 0, sizeof(zeros));
    for (int r = 0; r < N / PKE_COEFFS_PER_REG; r++)
      eccPKEWriteBuf(PKE_SLOT1 + r, zeros, PKE_COEFFS_PER_REG);
  }
  for (int j = 0; j < L; j++) {
    /* Expand y_j and NTT */
    poly_expand_mask(poly0, rho_prime, kappa + (uint16_t)j);
    ntt(poly0);
    /* Expand A[row][j] (already in NTT domain) */
    poly_rej_ntt(poly1, rho, (uint8_t)row, (uint8_t)j);
    /* w_hat += A_ij * y_hat_j (accumulated in PKE_SLOT1) */
    pointwise_acc_pke(poly1, poly0);
  }
}

/* ------------------------------------------------------------------ */
/*  ML-DSA-65 Signing (FIPS 204, Algorithm 7)                          */
/*                                                                      */
/*  Three-pass approach to minimize memory:                             */
/*    Pass 1: Compute w = Ay, hash w1 → c_tilde.                      */
/*    Pass 2: Compute z = y + c*s1, check ||z||∞ < γ1 − β.           */
/*    Pass 3: Recompute w, check r0 = LowBits(w − c*s2),              */
/*            compute c*t0, make hint.                                  */
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
  /* PKE_SLOT0 (regs 4-19):  challenge c (persistent across passes)
   * PKE_SLOT1 (regs 20-35): w_hat accumulator / temp storage
   * Total stack ≈ 2048 + 208 + 64 + 64 + 48 + misc ≈ 2.5 KB */

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
  for (;;) {
    reject = 0;

    /* ==== Pass 1: Compute c_tilde = H(mu || w1Encode(w)) ==== */
    shake256_init(&shake_ctx);
    shake_update(&shake_ctx, mu, 64);

    for (int i = 0; i < K; i++) {
      /* Compute w_hat[i] → PKE_SLOT1, using poly0/poly1 as scratch */
      compute_w_hat_row(sk_rho(sk), rho_prime, kappa, i, poly0, poly1);
      /* Load w_hat from PKE_SLOT1, INTT to get w[i] */
      pke_load_poly(poly0, PKE_SLOT1);
      invntt(poly0);
      poly_caddq(poly0);
      /* Encode w1 = HighBits(w[i]) directly into SHAKE */
      w1_encode_update(&shake_ctx, poly0);
    }

    shake_finalize(&shake_ctx);
    shake_squeeze(&shake_ctx, c_tilde, MLDSA_C_TILDE_BYTES);

    /* ==== Compute challenge c ==== */
    poly_challenge(poly0, c_tilde);
    poly_caddq(poly0);  /* ensure [0, Q) for Kronecker packing */
    /* Store c in PKE slot 0 (normal domain, regs 4-19) */
    pke_store_poly(PKE_SLOT0, poly0);

    /* ==== Pass 2: Compute z = y + c*s1, check bounds ==== */
    /* Also encode z into signature as we go. */
    for (int j = 0; j < L; j++) {
      /* y_j → poly0, save to PKE_SLOT1 */
      poly_expand_mask(poly0, rho_prime, kappa + (uint16_t)j);
      pke_store_poly(PKE_SLOT1, poly0);
      /* s1_j → poly0 */
      unpack_eta(poly0, sk_s1(sk, j));
      poly_caddq(poly0);  /* ensure [0, Q) */
      /* poly1 = c * s1_j (c read from PKE_SLOT0 inside) */
      poly_mul_challenge(poly1, poly0);
      /* Restore y_j from PKE_SLOT1 */
      pke_load_poly(poly0, PKE_SLOT1);
      /* z_j = y_j + c*s1_j */
      poly_add(poly0, poly0, poly1);
      poly_reduce(poly0);
      /* Center and check ||z_j||∞ >= gamma1 - beta → reject */
      for (int n = 0; n < N; n++) {
        if (poly0[n] > Q / 2) poly0[n] -= Q;
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
    if (reject) { kappa += L; continue; }

    /* ==== Pass 3: Recompute w, check r0 & ct0, make hint ==== */
    {
      uint8_t *hint_buf = sig + MLDSA_C_TILDE_BYTES
                        + L * MLDSA_POLYZ_PACKEDBYTES;
      int hint_count = 0;
      memset(hint_buf, 0, OMEGA + K);

      for (int i = 0; i < K; i++) {
        /* Recompute w[i] → PKE_SLOT1 */
        compute_w_hat_row(sk_rho(sk), rho_prime, kappa, i, poly0, poly1);
        pke_load_poly(poly0, PKE_SLOT1);
        invntt(poly0);
        poly_caddq(poly0);
        /* poly0 = w[i], save to PKE_SLOT1 for later */
        pke_store_poly(PKE_SLOT1, poly0);

        /* Compute c*s2[i]: s2 → poly0, result → poly1 */
        unpack_eta(poly0, sk_s2(sk, i));
        poly_caddq(poly0);  /* ensure [0, Q) */
        poly_mul_challenge(poly1, poly0);
        /* poly1 = c*s2[i] */

        /* Restore w[i], compute r = w[i] - c*s2[i] */
        pke_load_poly(poly0, PKE_SLOT1);
        poly_sub(poly0, poly0, poly1);
        poly_reduce(poly0);
        /* poly0 = r = w[i] - c*s2[i] */

        /* Check ||LowBits(r)||∞ >= gamma2 - beta → reject */
        for (int n = 0; n < N; n++) {
          int32_t r0 = low_bits(poly0[n]);
          if (r0 >= (int32_t)(GAMMA2 - BETA_B) ||
              r0 <= -(int32_t)(GAMMA2 - BETA_B)) {
            reject = 1;
            break;
          }
        }
        if (reject) break;

        /* Save r = w-cs2 to PKE_SLOT1, compute c*t0[i] */
        pke_store_poly(PKE_SLOT1, poly0);
        unpack_t0(poly0, sk_t0(sk, i));
        poly_caddq(poly0);  /* ensure [0, Q) */
        poly_mul_challenge(poly1, poly0);
        /* poly1 = c*t0[i] */
        /* Restore r = w - cs2 */
        pke_load_poly(poly0, PKE_SLOT1);

        /* Check ||c*t0||∞ < gamma2 */
        for (int n = 0; n < N; n++) {
          int32_t v = poly1[n];
          if (v > Q / 2) v -= Q;
          if (v >= (int32_t)GAMMA2 || v <= -(int32_t)GAMMA2) {
            reject = 1;
            break;
          }
        }
        if (reject) break;

        /* Make hint: h[i] = MakeHint(-ct0, w - cs2 + ct0) */
        /* poly0 = w - cs2, poly1 = ct0 */
        for (int n = 0; n < N; n++) {
          int32_t neg_ct0 = freeze(-poly1[n]);
          int32_t w_cs2_ct0 = freeze(poly0[n] + poly1[n]);
          if (make_hint(neg_ct0, w_cs2_ct0)) {
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
    /* c_tilde is at the beginning */
    memcpy(sig, c_tilde, MLDSA_C_TILDE_BYTES);
    /* z is already encoded in pass 2 */
    /* hint is already encoded in pass 3 */

    *sig_len = MLDSA_SIG_BYTES;
    return 0;
  }
}

/* ------------------------------------------------------------------ */
/*  ML-DSA-65 Key Generation (FIPS 204, Algorithm 6)                   */
/*                                                                      */
/*  pk = (rho || t1_packed)                                             */
/*  sk = (rho || K || tr || s1_packed || s2_packed || t0_packed)        */
/*                                                                      */
/*  Memory: reuses 3 stack polys + PKE buffer.                          */
/*  t = NTT^{-1}(A_hat * NTT(s1)) + s2, then Power2Round.             */
/* ------------------------------------------------------------------ */

int ml_dsa_65_keygen(uint8_t *pk, uint8_t *sk, const uint8_t *seed) {
  poly poly0, poly1, poly2;
  SHA3_CTX_T shake_ctx;
  uint8_t rho[32], rho_prime[64], K_seed[32];

  /* Step 1: (rho, rho', K) = H(seed) — squeeze 32+64+32 = 128 bytes */
  shake256_init(&shake_ctx);
  shake_update(&shake_ctx, seed, 32);
  {
    uint8_t dom[2] = { (uint8_t)K, (uint8_t)L };
    shake_update(&shake_ctx, dom, 2);
  }
  shake_finalize(&shake_ctx);
  shake_squeeze(&shake_ctx, rho, 32);
  shake_squeeze(&shake_ctx, rho_prime, 64);
  shake_squeeze(&shake_ctx, K_seed, 32);

  /* Copy rho and K into sk */
  memcpy(sk, rho, 32);           /* sk[0..31]  = rho */
  memcpy(sk + 32, K_seed, 32);   /* sk[32..63] = K   */

  /* Step 2: Generate s1, s2 using rho' and pack into sk */
  for (int j = 0; j < L; j++) {
    poly_rej_bounded(poly0, rho_prime, (uint16_t)j);
    pack_eta(sk + 128 + j * MLDSA_POLYETA_PACKEDBYTES, poly0);
  }
  for (int i = 0; i < K; i++) {
    poly_rej_bounded(poly0, rho_prime, (uint16_t)(L + i));
    pack_eta(sk + 128 + L * MLDSA_POLYETA_PACKEDBYTES
             + i * MLDSA_POLYETA_PACKEDBYTES, poly0);
  }

  /* Copy rho to pk */
  memcpy(pk, rho, 32);

  /* Step 3: Compute t = A*s1 + s2, row by row.
   * For each row i, compute t[i], then Power2Round → (t1, t0).
   * Pack t1 into pk, t0 into sk. */
  for (int i = 0; i < K; i++) {
    /* t_hat[i] = sum_j A_hat[i][j] * s1_hat[j] */
    memset(poly2, 0, N * sizeof(int32_t));
    for (int j = 0; j < L; j++) {
      unpack_eta(poly0, sk + 128 + j * MLDSA_POLYETA_PACKEDBYTES);
      ntt(poly0);
      poly_rej_ntt(poly1, rho, (uint8_t)i, (uint8_t)j);
      poly_pointwise_acc(poly2, poly1, poly0);
    }
    invntt(poly2);

    /* Add s2[i] */
    unpack_eta(poly0, sk + 128 + L * MLDSA_POLYETA_PACKEDBYTES
               + i * MLDSA_POLYETA_PACKEDBYTES);
    poly_add(poly2, poly2, poly0);
    poly_caddq(poly2);

    /* Power2Round: t = t1 * 2^D + t0 */
    for (int n = 0; n < N; n++) {
      int32_t t0_coeff;
      poly0[n] = power2round(&t0_coeff, poly2[n]); /* poly0 = t1 */
      poly1[n] = t0_coeff;                          /* poly1 = t0 */
    }

    /* Pack t1 into pk */
    pack_t1(pk + 32 + i * 320, poly0);
    /* Pack t0 into sk */
    pack_t0(sk + 128 + (L + K) * MLDSA_POLYETA_PACKEDBYTES
            + i * MLDSA_POLYT0_PACKEDBYTES, poly1);
  }

  /* Step 4: Compute tr = H(pk) and store in sk[64..127] */
  shake256_init(&shake_ctx);
  shake_update(&shake_ctx, pk, MLDSA_PK_BYTES);
  shake_finalize(&shake_ctx);
  shake_squeeze(&shake_ctx, sk + 64, MLDSA_TRBYTES);

  return 0;
}

/* ------------------------------------------------------------------ */
/*  Unpacking helpers for verification                                 */
/* ------------------------------------------------------------------ */

/* Unpack z polynomial (gamma1=2^19, 20 bits per coeff). */
static void unpack_z(int32_t a[N], const uint8_t *buf) {
  for (int i = 0; i < N / 2; i++) {
    uint32_t t0 = (uint32_t)buf[5*i+0] | ((uint32_t)buf[5*i+1] << 8) |
                  ((uint32_t)(buf[5*i+2] & 0x0F) << 16);
    uint32_t t1 = ((uint32_t)buf[5*i+2] >> 4) | ((uint32_t)buf[5*i+3] << 4) |
                  ((uint32_t)buf[5*i+4] << 12);
    a[2*i]   = (int32_t)(GAMMA1 - t0);
    a[2*i+1] = (int32_t)(GAMMA1 - t1);
  }
}

/* Unpack t1 polynomial (10 bits per coeff). */
static void unpack_t1(int32_t a[N], const uint8_t *buf) {
  for (int i = 0; i < N / 4; i++) {
    a[4*i+0] = (int32_t)(((uint32_t)buf[5*i+0] | ((uint32_t)buf[5*i+1] << 8)) & 0x3FF);
    a[4*i+1] = (int32_t)((((uint32_t)buf[5*i+1] >> 2) | ((uint32_t)buf[5*i+2] << 6)) & 0x3FF);
    a[4*i+2] = (int32_t)((((uint32_t)buf[5*i+2] >> 4) | ((uint32_t)buf[5*i+3] << 4)) & 0x3FF);
    a[4*i+3] = (int32_t)((((uint32_t)buf[5*i+3] >> 6) | ((uint32_t)buf[5*i+4] << 2)) & 0x3FF);
  }
}

/* UseHint: recover w1 from hint and r. FIPS 204, Algorithm 35. */
static int32_t use_hint(int32_t hint, int32_t r) {
  int32_t r0;
  int32_t r1 = decompose(&r0, r);
  if (hint == 0) return r1;
  if (r0 > 0) return (r1 + 1) & 15;
  return (r1 - 1) & 15;
}

/* ------------------------------------------------------------------ */
/*  ML-DSA-65 Verification (FIPS 204, Algorithm 3)                     */
/* ------------------------------------------------------------------ */

int ml_dsa_65_verify(const uint8_t *msg, size_t msg_len,
                     const uint8_t *sig, size_t sig_len,
                     const uint8_t *ctx, size_t ctx_len,
                     const uint8_t *pk) {
  poly poly0, poly1, poly2;
  SHA3_CTX_T shake_ctx;
  uint8_t mu[64], tr[MLDSA_TRBYTES];
  uint8_t c_tilde_check[MLDSA_C_TILDE_BYTES];

  if (sig_len != MLDSA_SIG_BYTES) return -1;
  if (ctx_len > 255) return -1;

  const uint8_t *rho = pk;
  const uint8_t *c_tilde = sig;
  const uint8_t *z_bytes = sig + MLDSA_C_TILDE_BYTES;
  const uint8_t *hint_buf = sig + MLDSA_C_TILDE_BYTES
                          + L * MLDSA_POLYZ_PACKEDBYTES;

  /* Decode and check hint format */
  {
    int prev = 0;
    for (int i = 0; i < K; i++) {
      int cur = (int)hint_buf[OMEGA + i];
      if (cur < prev || cur > OMEGA) return -2;
      prev = cur;
    }
  }

  /* Check ||z||_inf < gamma1 - beta */
  for (int j = 0; j < L; j++) {
    unpack_z(poly0, z_bytes + j * MLDSA_POLYZ_PACKEDBYTES);
    for (int n = 0; n < N; n++) {
      if (poly0[n] >= GAMMA1 - BETA_B || poly0[n] <= -(GAMMA1 - BETA_B))
        return -3;
    }
  }

  /* Compute tr = H(pk) */
  shake256_init(&shake_ctx);
  shake_update(&shake_ctx, pk, MLDSA_PK_BYTES);
  shake_finalize(&shake_ctx);
  shake_squeeze(&shake_ctx, tr, MLDSA_TRBYTES);

  /* Compute mu = H(tr || 0x00 || ctx_len || ctx || msg) */
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

  /* Compute c = SampleInBall(c_tilde) */
  poly_challenge(poly0, c_tilde);
  poly_caddq(poly0);  /* ensure [0, Q) for Kronecker packing */
  /* Store c in PKE slot 0 (normal domain) */
  pke_store_poly(PKE_SLOT0, poly0);

  /* Compute w'_approx and hash w1' */
  shake256_init(&shake_ctx);
  shake_update(&shake_ctx, mu, 64);

  for (int i = 0; i < K; i++) {
    /* w'_approx[i] = A[i]*z - c*(t1[i] << D) */

    /* Compute sum_j A_hat[i][j] * z_hat[j] in NTT domain, then INTT */
    memset(poly2, 0, N * sizeof(int32_t));
    for (int j = 0; j < L; j++) {
      unpack_z(poly0, z_bytes + j * MLDSA_POLYZ_PACKEDBYTES);
      ntt(poly0);
      poly_rej_ntt(poly1, rho, (uint8_t)i, (uint8_t)j);
      poly_pointwise_acc(poly2, poly1, poly0);
    }
    invntt(poly2);
    poly_caddq(poly2);
    /* poly2 = A[i]*z in normal domain */

    /* Compute c * (t1[i] << D) via Kronecker (c is in PKE_SLOT0) */
    unpack_t1(poly0, pk + 32 + i * MLDSA_POLYT1_PACKEDBYTES);
    for (int n = 0; n < N; n++)
      poly0[n] <<= D_BITS;
    poly_reduce(poly0);
    poly_caddq(poly0);  /* ensure [0, Q) for Kronecker packing */
    poly_mul_challenge(poly1, poly0);
    /* poly1 = c * (t1[i] << D) */

    /* w'_approx = Az - c*(t1<<D) */
    poly_sub(poly2, poly2, poly1);
    poly_reduce(poly2);
    poly_caddq(poly2);

    /* Apply hint to get w1' */
    {
      int h_start = (i == 0) ? 0 : (int)hint_buf[OMEGA + i - 1];
      int h_end   = (int)hint_buf[OMEGA + i];
      memset(poly0, 0, N * sizeof(int32_t));
      for (int h = h_start; h < h_end; h++)
        poly0[hint_buf[h]] = 1;
      for (int n = 0; n < N; n++)
        poly1[n] = use_hint(poly0[n], poly2[n]);
    }
    w1_pack_update(&shake_ctx, poly1);
  }

  shake_finalize(&shake_ctx);
  shake_squeeze(&shake_ctx, c_tilde_check, MLDSA_C_TILDE_BYTES);

  /* Compare c_tilde */
  if (memcmp(c_tilde, c_tilde_check, MLDSA_C_TILDE_BYTES) != 0)
    return -4;

  return 0;
}

/* ------------------------------------------------------------------ */
/*  Self-test (FIPS 140-3 style KAT + round-trip)                      */
/* ------------------------------------------------------------------ */

int ml_dsa_65_selftest(void) {
  /* ACVP tcId 26 seed */
  static const uint8_t seed[32] = {
    0x1B, 0xD6, 0x7D, 0xC7, 0x82, 0xB2, 0x95, 0x8E,
    0x18, 0x9E, 0x31, 0x5C, 0x04, 0x0D, 0xD1, 0xF6,
    0x4C, 0x8A, 0xB2, 0x32, 0xA6, 0xA1, 0x70, 0xE1,
    0xA7, 0xA5, 0x2C, 0x33, 0xF1, 0x08, 0x51, 0xB1
  };

  /* Expected SHA3-256 digests of pk and sk */
  static const uint8_t expected_pk_hash[32] = {
    0xaa, 0xa0, 0x7f, 0x58, 0x6d, 0x78, 0xb6, 0x7b,
    0x96, 0x4d, 0xe8, 0xde, 0xf0, 0xdf, 0x7f, 0x34,
    0xa6, 0xc1, 0x60, 0xf1, 0x10, 0xba, 0x70, 0x1a,
    0x7c, 0x1a, 0x28, 0xb9, 0xba, 0x2f, 0x8b, 0xa6
  };
  static const uint8_t expected_sk_hash[32] = {
    0x0a, 0xd8, 0xc5, 0x37, 0x1e, 0x61, 0xcd, 0x02,
    0x6e, 0x0d, 0xad, 0x72, 0xdf, 0xc3, 0x74, 0x08,
    0x40, 0x18, 0x7e, 0x20, 0x94, 0xc2, 0x7f, 0x91,
    0x5e, 0xb5, 0xce, 0xf8, 0xfd, 0x19, 0x12, 0x6e
  };

  /* Expected SHA3-256 digest of signature (deterministic, msg below, no ctx) */
  static const uint8_t expected_sig_hash[32] = {
    0xf9, 0xad, 0x98, 0xe1, 0x6f, 0xdd, 0x68, 0x0c,
    0xa7, 0x85, 0x86, 0x51, 0xb0, 0xef, 0x3f, 0x16,
    0x34, 0x36, 0xa0, 0x74, 0xdf, 0xf9, 0x6e, 0x63,
    0x2f, 0xb7, 0x9b, 0x13, 0x64, 0x50, 0xfb, 0x31
  };

  static uint8_t pk[MLDSA_PK_BYTES];
  static uint8_t sk[MLDSA_SK_BYTES];
  static uint8_t sig[MLDSA_SIG_BYTES];
  uint8_t digest[32];
  size_t sig_len = 0;

  static const uint8_t msg[] = "test message for ML-DSA-65";
  static const size_t msg_len = sizeof(msg) - 1;

  /* 1. KeyGen KAT */
  if (ml_dsa_65_keygen(pk, sk, seed) != 0)
    return -1;

  sha3_256_raw(pk, MLDSA_PK_BYTES, digest);
  if (memcmp(digest, expected_pk_hash, 32) != 0)
    return -2;

  sha3_256_raw(sk, MLDSA_SK_BYTES, digest);
  if (memcmp(digest, expected_sk_hash, 32) != 0)
    return -3;

  /* 2. Sign KAT (deterministic) */
  if (ml_dsa_65_sign(sig, &sig_len, msg, msg_len, NULL, 0, sk) != 0)
    return -4;

  if (sig_len != MLDSA_SIG_BYTES)
    return -5;

  sha3_256_raw(sig, sig_len, digest);
  if (memcmp(digest, expected_sig_hash, 32) != 0)
    return -6;

  /* 3. Verify round-trip */
  if (ml_dsa_65_verify(msg, msg_len, sig, sig_len, NULL, 0, pk) != 0)
    return -7;

  /* 4. Verify rejects wrong message */
  {
    static const uint8_t bad[] = "wrong message";
    if (ml_dsa_65_verify(bad, sizeof(bad) - 1, sig, sig_len, NULL, 0, pk) == 0)
      return -8;
  }

  return 0;
}

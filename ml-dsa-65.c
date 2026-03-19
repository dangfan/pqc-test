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
/*  Hardware MMM setup for Q                                           */
/* ------------------------------------------------------------------ */

static void pke_mmm_setup(void) {
  uint32_t mod = Q;
  rsaPKESetLen(1);
  eccPKEWriteBuf(0, &mod, 1);            /* reg 0 = Q */
  uint32_t mc[2];
  eccCalMc64(mc, &mod);
  rsaPKEWriteMc(mc);
}

/* Hardware modular multiply: a * b * R^{-1} mod Q */
static int32_t hw_mon_mul(int32_t a, int32_t b) {
  uint32_t ua = (uint32_t)freeze(a);
  uint32_t ub = (uint32_t)freeze(b);
  uint32_t result;
  eccPKEWriteBuf(1, &ua, 1);
  eccPKEWriteBuf(2, &ub, 1);
  sm2MonMul(3, 1, 2);
  eccPKEReadBuf(&result, 3, 1);
  return (int32_t)result;
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
/*  Uses poly0 for y_hat_j, poly1 for A_hat_ij, result in w_hat.      */
/*  w_hat is accumulated in PKE_SLOT1 to free stack poly buffers.      */
/* ------------------------------------------------------------------ */

static void compute_w_hat_row(int32_t w_hat[N],
                              const uint8_t *rho,
                              const uint8_t *rho_prime,
                              uint16_t kappa, int row,
                              int32_t poly0[N], int32_t poly1[N]) {
  memset(w_hat, 0, N * sizeof(int32_t));
  for (int j = 0; j < L; j++) {
    /* Expand y_j and NTT */
    poly_expand_mask(poly0, rho_prime, kappa + (uint16_t)j);
    ntt(poly0);
    /* Expand A[row][j] (already in NTT domain) */
    poly_rej_ntt(poly1, rho, (uint8_t)row, (uint8_t)j);
    /* w_hat += A_ij * y_hat_j */
    poly_pointwise_acc(w_hat, poly1, poly0);
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
  poly poly0, poly1, poly2;             /* 3 × 1024 = 3072 bytes */
  SHA3_CTX_T shake_ctx;                 /* ~208 bytes */
  uint8_t mu[64];
  uint8_t rho_prime[64];
  uint8_t c_tilde[MLDSA_C_TILDE_BYTES]; /* 48 bytes */
  uint8_t w1_packed[MLDSA_POLYW1_PACKEDBYTES]; /* 128 bytes */
  uint16_t kappa;
  int reject;
  /* Total stack ≈ 3072 + 208 + 64 + 64 + 48 + 128 + misc ≈ 3.6 KB */

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
      /* Compute w_hat[i] in poly2, using poly0/poly1 as scratch */
      compute_w_hat_row(poly2, sk_rho(sk), rho_prime, kappa, i,
                        poly0, poly1);
      /* INTT to get w[i] in normal domain */
      invntt(poly2);
      poly_caddq(poly2);
      /* Extract w1 = HighBits(w[i]) and encode */
      for (int n = 0; n < N; n++)
        poly0[n] = high_bits(poly2[n]);
      pack_w1(w1_packed, poly0);
      shake_update(&shake_ctx, w1_packed, MLDSA_POLYW1_PACKEDBYTES);
    }

    shake_finalize(&shake_ctx);
    shake_squeeze(&shake_ctx, c_tilde, MLDSA_C_TILDE_BYTES);

    /* ==== Compute challenge c and NTT(c) ==== */
    poly_challenge(poly0, c_tilde);
    ntt(poly0);
    /* Store c_hat in PKE slot 0 (regs 4-19) */
    pke_store_poly(PKE_SLOT0, poly0);

    /* ==== Pass 2: Compute z = y + c*s1, check bounds ==== */
    /* Also encode z into signature as we go. */
    for (int j = 0; j < L; j++) {
      /* y_j */
      poly_expand_mask(poly0, rho_prime, kappa + (uint16_t)j);
      /* c_hat * NTT(s1_j) */
      unpack_eta(poly1, sk_s1(sk, j));
      ntt(poly1);
      pke_load_poly(poly2, PKE_SLOT0); /* c_hat */
      poly_pointwise(poly2, poly2, poly1);
      invntt(poly2);
      /* z_j = y_j + INTT(c_hat * s1_hat_j) */
      poly_add(poly0, poly0, poly2);
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
        /* Recompute w[i] */
        compute_w_hat_row(poly2, sk_rho(sk), rho_prime, kappa, i,
                          poly0, poly1);
        invntt(poly2);
        poly_caddq(poly2);
        /* poly2 = w[i] */

        /* Compute c*s2[i] */
        unpack_eta(poly0, sk_s2(sk, i));
        ntt(poly0);
        pke_load_poly(poly1, PKE_SLOT0); /* c_hat */
        poly_pointwise(poly0, poly1, poly0);
        invntt(poly0);
        /* poly0 = c*s2[i] */

        /* r = w[i] - c*s2[i] */
        poly_sub(poly1, poly2, poly0);
        poly_reduce(poly1);
        /* Check ||LowBits(r)||∞ >= gamma2 - beta → reject */
        for (int n = 0; n < N; n++) {
          int32_t r0 = low_bits(poly1[n]);
          if (r0 >= (int32_t)(GAMMA2 - BETA_B) ||
              r0 <= -(int32_t)(GAMMA2 - BETA_B)) {
            reject = 1;
            break;
          }
        }
        if (reject) break;

        /* Compute c*t0[i] — save poly1 (w-cs2) to PKE slot 1 temporarily */
        pke_store_poly(PKE_SLOT1, poly1);
        unpack_t0(poly0, sk_t0(sk, i));
        ntt(poly0);
        pke_load_poly(poly1, PKE_SLOT0); /* c_hat */
        poly_pointwise(poly0, poly1, poly0);
        invntt(poly0);
        poly_reduce(poly0);
        /* Restore poly1 = w - cs2 */
        pke_load_poly(poly1, PKE_SLOT1);
        /* Check ||c*t0||∞ < gamma2 */
        for (int n = 0; n < N; n++) {
          int32_t v = poly0[n];
          if (v > Q / 2) v -= Q;
          if (v >= (int32_t)GAMMA2 || v <= -(int32_t)GAMMA2) {
            reject = 1;
            break;
          }
        }
        if (reject) break;

        /* Make hint: h[i] = MakeHint(-ct0, w - cs2 + ct0) */
        /* poly1 = w - cs2 (already computed), poly0 = ct0 */
        for (int n = 0; n < N; n++) {
          int32_t neg_ct0 = freeze(-poly0[n]);
          int32_t w_cs2_ct0 = freeze(poly1[n] + poly0[n]);
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
  uint8_t w1_packed[MLDSA_POLYW1_PACKEDBYTES];

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

  /* Compute c = SampleInBall(c_tilde) and NTT(c) */
  poly_challenge(poly0, c_tilde);
  ntt(poly0);
  /* Store c_hat in PKE slot 0 */
  pke_store_poly(PKE_SLOT0, poly0);

  /* Compute w'_approx and hash w1' */
  shake256_init(&shake_ctx);
  shake_update(&shake_ctx, mu, 64);

  for (int i = 0; i < K; i++) {
    /* w'_approx[i] = NTT^{-1}(A_hat[i] * z_hat - c_hat * NTT(t1[i] << D)) */

    /* Compute sum_j A_hat[i][j] * z_hat[j] */
    memset(poly2, 0, N * sizeof(int32_t));
    for (int j = 0; j < L; j++) {
      unpack_z(poly0, z_bytes + j * MLDSA_POLYZ_PACKEDBYTES);
      ntt(poly0);
      poly_rej_ntt(poly1, rho, (uint8_t)i, (uint8_t)j);
      poly_pointwise_acc(poly2, poly1, poly0);
    }

    /* Compute c_hat * NTT(t1[i] << D) */
    unpack_t1(poly0, pk + 32 + i * MLDSA_POLYT1_PACKEDBYTES);
    for (int n = 0; n < N; n++)
      poly0[n] <<= D_BITS;
    ntt(poly0);
    pke_load_poly(poly1, PKE_SLOT0); /* c_hat */
    poly_pointwise(poly0, poly1, poly0);

    /* w'_approx = Az_hat - ct1_hat */
    poly_sub(poly2, poly2, poly0);
    invntt(poly2);
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
    pack_w1(w1_packed, poly1);
    shake_update(&shake_ctx, w1_packed, MLDSA_POLYW1_PACKEDBYTES);
  }

  shake_finalize(&shake_ctx);
  shake_squeeze(&shake_ctx, c_tilde_check, MLDSA_C_TILDE_BYTES);

  /* Compare c_tilde */
  if (memcmp(c_tilde, c_tilde_check, MLDSA_C_TILDE_BYTES) != 0)
    return -4;

  return 0;
}

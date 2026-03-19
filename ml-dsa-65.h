#ifndef _ML_DSA_65_H_
#define _ML_DSA_65_H_

#include <stddef.h>
#include <stdint.h>

/* ML-DSA-65 (FIPS 204) parameters */
#define MLDSA_Q          8380417
#define MLDSA_N          256
#define MLDSA_D          13
#define MLDSA_TAU        49
#define MLDSA_LAMBDA     192
#define MLDSA_GAMMA1     (1 << 19)
#define MLDSA_GAMMA2     ((MLDSA_Q - 1) / 32)  /* 261888 */
#define MLDSA_K          6
#define MLDSA_L          5
#define MLDSA_ETA        4
#define MLDSA_BETA       (MLDSA_TAU * MLDSA_ETA) /* 196 */
#define MLDSA_OMEGA      55

/* Derived sizes (bytes) */
#define MLDSA_SEEDBYTES       32
#define MLDSA_CRHBYTES        64
#define MLDSA_TRBYTES         64
#define MLDSA_C_TILDE_BYTES   48   /* lambda / 4 */
#define MLDSA_POLYW1_PACKEDBYTES  128  /* 256 * 4 / 8 */

/* Secret key component packed sizes */
#define MLDSA_POLYETA_PACKEDBYTES   128  /* eta=4: 4 bits per coeff, 256*4/8 */
#define MLDSA_POLYT0_PACKEDBYTES    416  /* 13 bits per coeff, 256*13/8 */
#define MLDSA_POLYZ_PACKEDBYTES     640  /* gamma1=2^19: 20 bits, 256*20/8 */
#define MLDSA_POLYT1_PACKEDBYTES    320  /* 10 bits per coeff, 256*10/8 */

/* Key sizes */
#define MLDSA_PK_BYTES  (MLDSA_SEEDBYTES + MLDSA_K * MLDSA_POLYT1_PACKEDBYTES)
/* 32 + 6*320 = 1952 bytes */

#define MLDSA_SK_BYTES (MLDSA_SEEDBYTES + MLDSA_SEEDBYTES + MLDSA_TRBYTES \
                        + MLDSA_L * MLDSA_POLYETA_PACKEDBYTES             \
                        + MLDSA_K * MLDSA_POLYETA_PACKEDBYTES             \
                        + MLDSA_K * MLDSA_POLYT0_PACKEDBYTES)
/* 32 + 32 + 64 + 5*128 + 6*128 + 6*416 = 4000 bytes */

/* Signature size */
#define MLDSA_SIG_BYTES (MLDSA_C_TILDE_BYTES                   \
                         + MLDSA_L * MLDSA_POLYZ_PACKEDBYTES    \
                         + MLDSA_OMEGA + MLDSA_K)
/* 48 + 5*640 + 55 + 6 = 3309 bytes */

/**
 * ML-DSA-65 signing (FIPS 204, Algorithm 7 – ML-DSA.Sign).
 *
 * @param sig     Output signature buffer, MLDSA_SIG_BYTES bytes.
 * @param sig_len Set to the actual signature length on success.
 * @param msg     Message to sign.
 * @param msg_len Length of message.
 * @param ctx     Context string (may be NULL if ctx_len == 0).
 * @param ctx_len Context string length (0..255).
 * @param sk      Secret key, MLDSA_SK_BYTES bytes.
 *
 * @return 0 on success, negative on failure.
 */
int ml_dsa_65_sign(uint8_t *sig, size_t *sig_len,
                   const uint8_t *msg, size_t msg_len,
                   const uint8_t *ctx, size_t ctx_len,
                   const uint8_t *sk);

/**
 * ML-DSA-65 key generation (FIPS 204, Algorithm 6 – ML-DSA.KeyGen).
 *
 * @param pk      Output public key buffer.
 * @param sk      Output secret key buffer, MLDSA_SK_BYTES bytes.
 * @param seed    32-byte random seed (xi).
 *
 * @return 0 on success, negative on failure.
 */
int ml_dsa_65_keygen(uint8_t *pk, uint8_t *sk, const uint8_t *seed);

/**
 * ML-DSA-65 verification (FIPS 204, Algorithm 3 – ML-DSA.Verify).
 *
 * @param msg     Message that was signed.
 * @param msg_len Length of message.
 * @param sig     Signature, MLDSA_SIG_BYTES bytes.
 * @param sig_len Signature length (must equal MLDSA_SIG_BYTES).
 * @param ctx     Context string (may be NULL if ctx_len == 0).
 * @param ctx_len Context string length (0..255).
 * @param pk      Public key, MLDSA_PK_BYTES bytes.
 *
 * @return 0 on success (valid signature), negative on failure.
 */
int ml_dsa_65_verify(const uint8_t *msg, size_t msg_len,
                     const uint8_t *sig, size_t sig_len,
                     const uint8_t *ctx, size_t ctx_len,
                     const uint8_t *pk);

/* Debug: test NTT/INTT round-trip. Returns 0 on success. */
int ml_dsa_65_selftest(void);

#endif /* _ML_DSA_65_H_ */

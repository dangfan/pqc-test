#ifndef _SE_H_
#define _SE_H_

#include <stdint.h>

#define PRIME 0

// external functions
////////////////////////////////////////////////////////////////////////////////////////////////
/// PKE Register
/// 64 bytes * 48, big endian
/// A larger (than 32-byte) number can cross registers, e.g., 2^1024 - 1 can be stored as
/// 0: FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
///    FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
/// 1: FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
///    FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
///
/// NOTE: The modulus should be always stored in reg 0
////////////////////////////////////////////////////////////////////////////////////////////////

/**
 * Set the length of big numbers in the PKE co-processor.
 *
 * @param len Length in 4-byte, e.g., 64 for 256 bytes (2048 bit). Should be no greater than 64.
 *
 * @return 0 for success.
 */
int rsaPKESetLen(int len);

/**
 * Compute the 64-bit montgomery constant.
 *
 * @param mc Target data buffer to store the result. Should be uint32_t[2].
 * @param p  The pointer to the start of last eight bytes of the modulus.
 *
 * @return 0 for success.
 */
int eccCalMc64(uint32_t *mc, const uint32_t *p);

/**
 * Set the 64-bit montgomery constant.
 *
 * @param buf The constant to be set.
 *
 * @return 0 for success.
 */
int rsaPKEWriteMc(const uint32_t *buf);

/**
 * Write data to a PKE register.
 *
 * @param idx   Index of the target register.
 * @param buf   Data to write.
 * @param len   Length of data in 4-bytes.
 *
 * @return 0 for success.
 */
int eccPKEWriteBuf(int idx, const uint32_t *buf, int len);

/**
 * Read data from a PKE register.
 *
 * @param buf   Data to store the value.
 * @param idx   Index of the source register.
 * @param len   Length of data in 4-bytes.
 *
 * @return 0 for success.
 */
int eccPKEReadBuf(uint32_t *buf, int idx, int len);

/**
 * Montgomery modular addition. Operators can be the same register.
 * monAdd(A, B) = A + B mod N.
 * N should be previously set.
 * N: reg 0
 *
 * @param r  Target register.
 * @param a  Operator 1.
 * @param b  Operator 2.
 *
 * @return 0 for success.
 */
int sm2MonAdd(int r, int a, int b);

/**
 * Montgomery modular subtraction. Operators can be the same register.
 * monSub(A, B) = A - B mod N.
 * N should be previously set.
 * N: reg 0
 *
 * @param r  Target register.
 * @param a  Operator 1.
 * @param b  Operator 2.
 *
 * @return 0 for success.
 */
int sm2MonSub(int r, int a, int b);

/**
 * Montgomery modular multiplication. Operators can be the same register.
 * monMul(A, B) = A * B / R mod N.
 * R & N should be previously set.
 * R: rsaPKEWriteMc
 * N: reg 0
 *
 * @param r  Target register.
 * @param a  Operator 1.
 * @param b  Operator 2.
 *
 * @return 0 for success.
 */
int sm2MonMul(int r, int a, int b);

/**
 * Self tests for cryptographic algorithms that use SE
 */
void self_test(void);

// extra functions
void se_start(void);
void se_stop(void);

#endif // _SE_H_

// -----------------------------------------------------------------------
// Collatz step function calculator 
// -----------------------------------------------------------------------

#include <stdio.h>
#include <stdint.h>
#ifdef __x86_64__
#include <x86intrin.h>
#endif
#include "collatz.h"

#define USE_TABLES 1		// use small precomputed tables for small numbers
#define STRIP_MSBITS 0		// when 1 use smaller numbers to prevent "don't care" overflows

#define THRESHOLD64  20		// 32 * ln(2)/ln(3)
#define THRESHOLD128 40		// 64 * ln(2)/ln(3)

typedef unsigned __int128 uint128_t;

// hi:lo = a * b
static inline uint64_t my_mul64(uint64_t a, uint64_t b, uint64_t * hi)
{
#ifdef __x86_64__
#if defined(__BMI2__)
#if 0
	// slower in the context used by collatz
	// mulx : 4 cycles latency
	// mul  : 3 cycles latency
        return _mulx_u64(a, b, (unsigned long long *) hi);
#else
	uint64_t lo;
 asm(" mulxq %[a],%[lo],%[hi]\n": [lo] "=r"(lo), [hi] "=r"(*hi): [a] "r"(a), "d"(b):);
	return lo;
#endif
#else
	uint64_t lo;
 asm(" mulq %[b]\n": "=a"(lo), "=d"(*hi): [a] "a"(a),[b] "rm"(b):"flags");
	return lo;
#endif
#else
#error "ARM processors : use mull"
#endif
}

// c::out = a + b + c
static inline uint8_t my_add64(uint8_t c, uint64_t a, uint64_t b, uint64_t * out)
{
#ifdef __x86_64__
	return _addcarry_u64(c, a, b, (unsigned long long *)out);
#else
#error "ARM processors : use adc"
#endif
}

// (a::b) >>= c
static inline void my_shr64(uint64_t * a, uint64_t * b, uint64_t c)
{
 asm(" shrdq %%cl, %1, %0\n shrq %%cl, %1\n": "=r"(*b), "=r"(*a): "0"(*b), "1"(*a), "c"(c):"flags");
}

static inline uint64_t bzh64(uint64_t n, uint64_t pos)
{
#ifdef __BMI2__
	return _bzhi_u64(n, pos);
#else
	return n & ((1ull << pos) - 1);
#endif
}

static inline uint128_t bzh128(uint128_t n, uint64_t pos)
{
	if (pos >= 64) {
		uint128_t t = bzh64((uint64_t) (n >> 64), pos - 64);
		return (n & (uint64_t) - 1) | (t << 64);
	}
	return bzh64((uint64_t) n, pos);
}

static inline uint64_t my_ctz(uint64_t n)
{
#if defined(__BMI1__)
	uint64_t t;
 asm(" tzcntq %1, %0\n": "=r"(t): "r"(n):"flags");
	return t;
#else
	if (n)
		return __builtin_ctzll(n);
	return 64;
#endif
}

static inline uint64_t my_clz(uint64_t n)
{
#ifdef __BMI1__
	uint64_t t;
 asm(" lzcntq %1, %0\n": "=r"(t): "r"(n):"flags");
	return t;
#else
	if (n)
		return __builtin_clzll(n);
	return 64;
#endif
}

static uint128_t mpz_get_ui128(mpz_srcptr gmp_z) __GMP_NOTHROW
{
	mp_ptr gmp_p = gmp_z->_mp_d;
	mp_size_t gmp_n = gmp_z->_mp_size;

	unsigned i = 0;
	unsigned bits = 0;
	uint128_t t = 0;
	while (i < gmp_n && bits < 128) {
		uint128_t s = gmp_p[i];
		t |= s << bits;
		bits += GMP_NUMB_BITS;
		i += 1;
	}
	return t;
}

static void mpz_set_ui128(mpz_t t, uint128_t v)
{
	if (v >> 64) {
		mpz_set_ui(t, (uint64_t) (v >> 64));
		mpz_mul_2exp(t, t, 64);
		mpz_add_ui(t, t, (uint64_t) v);
	} else {
		mpz_set_ui(t, (uint64_t) v);
	}
}

static inline uint64_t clz128(uint128_t t)
{
	uint64_t hi = (uint64_t) (t >> 64);
	if (hi) {
		return my_clz(hi);
	}
	uint64_t lo = (uint64_t) t;
	return 64 + my_clz(lo);
}

static inline uint64_t ctz128(uint128_t t)
{
	uint64_t lo = (uint64_t) t;
	if (lo) {
		return my_ctz(lo);
	}
	uint64_t hi = (uint64_t) (t >> 64);
	return 64 + my_ctz(hi);
}

static inline uint64_t clz64(uint64_t t)
{
	return my_clz(t);
}

static inline uint64_t ctz64(uint64_t t)
{
	return my_ctz(t);
}

static uint64_t helper(mpz_t, mpz_t, mpz_t, uint64_t);
static uint64_t helper64(uint64_t, uint64_t *, uint64_t *, uint64_t);
static uint64_t helper128(uint128_t, uint64_t *, uint64_t *, uint64_t k);
static uint64_t collatz64(uint64_t);
static uint64_t collatz128(uint128_t);

uint64_t fastest_collatz(mpz_t n)
{
	uint64_t L2, count = 0;
	mpz_t r, d, cc;
	mpz_init(r);
	mpz_init(d);
	mpz_init(cc);

	L2 = mpz_sizeinbase(n, 2) >> 1;
	while (L2 >= THRESHOLD128) {
		// gmp_printf("Test helper n=%Zd.\n", n);
		mpz_fdiv_r_2exp(r, n, L2);
		count += helper(r, d, cc, L2);
		mpz_mul(n, n, cc);
		mpz_add(n, n, d);
		mpz_fdiv_q_2exp(n, n, L2);
		L2 = mpz_sizeinbase(n, 2) >> 1;
	}

	mpz_clear(cc);
	mpz_clear(d);
	mpz_clear(r);

	count += collatz128(mpz_get_ui128(n));

	return count;
}

// This internal entrypoint
// can process numbers up to ~40 bits without overflow
static uint64_t collatz128(uint128_t n)
{
	uint64_t L2, count = 0;
	uint64_t d, cc;
	uint128_t r, mask;

	L2 = (128 - clz128(n)) >> 1;
	while (L2 >= THRESHOLD64) {
#if STRIP_MSBITS
		r = bzh128(n, L2);
#else
		r = n;
#endif
		count += helper128(r, &d, &cc, L2);
		if (clz128(n) + clz64(cc) < 128) {
			mpz_t tn, t;

			/* n * cc + d triggers a 128 bit overflow, use gmp */
			mpz_init(tn);
			mpz_init(t);
			mpz_set_ui128(tn, n);
			mpz_set_ui128(t, cc);
			mpz_mul(tn, tn, t);
			mpz_set_ui128(t, d);
			mpz_add(tn, tn, t);
			mpz_clear(t);
			mpz_fdiv_q_2exp(tn, tn, L2);
			if (mpz_sizeinbase(tn, 2) >= 128) {
				/* number is too large, continue with gmp */
				count += fastest_collatz(tn);
				mpz_clear(tn);
				return count;
			} else {
				n = mpz_get_ui128(tn);
				mpz_clear(tn);
			}
		} else {
			n *= cc;
			n += d;
			n >>= L2;
		}
		L2 = (128 - clz128(n)) >> 1;
	}

	count += collatz64((uint64_t) n);

	return count;
}

// This internal entrypoint
// can process numbers up to ~20 bits without overflow
static uint64_t collatz64(uint64_t n)
{
	uint64_t c, L2, count = 0;
	uint64_t hi, r, d, cc;

	while (n > 7) {
		L2 = (64 - clz64(n)) >> 1;
#if STRIP_MSBITS
		r = bzh64(n, L2);
#else
		r = n;
#endif
		count += helper64(r, &d, &cc, L2);
		// n = (n * cc + d) >> L2;
		n = my_mul64(n, cc, &hi);
		hi += my_add64(0, n, d, &n);
		my_shr64(&hi, &n, L2);
		if (hi) {
			/* number is too large, continue with uint128_t */
			uint128_t t = hi;
			t <<= 64;
			t |= n;
			return count + collatz128(t);
		}
	}

	static uint8_t table_a[8] = { 0, 0, 1, 7, 2, 5, 8, 16 };
	count += table_a[n];

	return count;
}

// n---> (3^c*n+d)/2^k
static uint64_t helper64(uint64_t n, uint64_t * d, uint64_t * cc, uint64_t k)
{

	uint64_t ans;
	uint64_t d1, d2, r, cc1, cc2;

	if (k == 2) {
		// harcoded small values (branch-free)
		// and avoids the recursivity overheads
#ifdef USE_TABLES
		static uint8_t prec_c[4] = { 1, 3, 3, 9 };
		static uint8_t prec_d[4] = { 0, 1, 2, 5 };
		static uint8_t prec_a[4] = { 2, 3, 3, 4 };

		n &= 3;
		*cc = prec_c[n];
		*d = prec_d[n];
		return prec_a[n];
#else
		// k1 = 1
		r = n & 1;
		cc1 = 1 + 2 * r;
		n = (n * cc1 + r) >> 1;
		// k2 = 1
		n &= 1;
		cc2 = 1 + 2 * n;
		*cc = cc1 * cc2;
		*d = (r * cc2) + (n << 1);
		return 2 + r + n;
#endif
	}

	if (k == 3) {
		// harcoded small values (branch-free)
		// and avoids the recursivity overheads
#ifdef USE_TABLES
		static uint8_t prec_c[8] = { 1, 9, 3, 9, 3, 3, 9, 27 };
		static uint8_t prec_d[8] = { 0, 7, 2, 5, 4, 1, 10, 19 };
		static uint8_t prec_a[8] = { 3, 5, 4, 5, 4, 4, 5, 6 };

		n &= 7;
		*cc = prec_c[n];
		*d = prec_d[n];
		return prec_a[n];
#else

		uint64_t cc0, d0;
		// k1 = 1
		r = n & 1;
		cc0 = 1 + 2 * r;
		d0 = r;
		ans = 1 + r;
		n = (n * cc0 + r) >> 1;
		// k2 = 2
		r = n & 1;
		cc1 = 1 + 2 * r;
		n = (n * cc1 + r) >> 1;
		n &= 1;
		cc2 = 1 + 2 * n;
		cc1 = cc1 * cc2;
		d1 = (r * cc2) + (n << 1);
		ans += 2 + r + n;
		*d = (d0 * cc1) + (d1 << 1);
		*cc = cc1 * cc0;
		return ans;
#endif
	}

#ifdef USE_TABLES
	if (k == 4) {
		static uint8_t prec_c[16] = { 1, 9, 9, 9, 3, 3, 9, 27, 3, 27, 3, 27, 9, 9, 27, 81 };
		static uint8_t prec_d[16] = { 0, 7, 14, 5, 4, 1, 10, 19, 8, 29, 2, 23, 20, 11, 38, 65 };
		static uint8_t prec_a[16] = { 4, 6, 6, 6, 5, 5, 6, 7, 5, 7, 5, 7, 6, 6, 7, 8 };
		n &= 15;
		*cc = prec_c[n];
		*d = prec_d[n];
		return prec_a[n];
	}

	if (k == 5) {
		static uint8_t prec_c[32] =
		    { 1, 27, 9, 9, 9, 9, 9, 81, 3, 81, 3, 27, 9, 9, 27, 81, 3, 9, 27, 27, 3, 3, 27, 27, 9, 27, 9, 81, 27, 27, 81,
			243
		};
		static uint8_t prec_d[32] =
		    { 0, 37, 14, 5, 28, 19, 10, 73, 8, 103, 2, 23, 20, 11, 38, 65, 16, 7, 58, 31, 4, 1, 46, 19, 40, 29, 22, 85, 76,
			49, 130, 211
		};
		static uint8_t prec_a[32] =
		    { 5, 8, 7, 7, 7, 7, 7, 9, 6, 9, 6, 8, 7, 7, 8, 9, 6, 7, 8, 8, 6, 6, 8, 8, 7, 8, 7, 9, 8, 8, 9, 10 };
		n &= 31;
		*cc = prec_c[n];
		*d = prec_d[n];
		return prec_a[n];
	}

	if (k == 6) {
static uint16_t prec_d[64] = {0, 37, 74, 47, 28, 19, 10, 73, 56, 103, 38, 23, 20, 11, 146, 65, 16, 53, 206, 125, 4, 1, 46, 19, 40, 29, 22, 287, 76, 49, 130, 211, 32, 143, 14, 5, 116, 89, 62, 251, 8, 341, 2, 101, 92, 65, 38, 227, 80, 7, 58, 31, 44, 35, 170, 89, 152, 119, 98, 85, 260, 179, 422, 665 }; 
static uint16_t prec_c[64] = {1, 27, 27, 27, 9, 9, 9, 81, 9, 81, 9, 27, 9, 9, 81, 81, 3, 27, 81, 81, 3, 3, 27, 27, 9, 27, 9, 243, 27, 27, 81, 243, 3, 81, 9, 9, 27, 27, 27, 243, 3, 243, 3, 81, 27, 27, 27, 243, 9, 9, 27, 27, 9, 9, 81, 81, 27, 81, 27, 81, 81, 81, 243, 729 }; 
static uint8_t prec_a[64] = {6, 9, 9, 9, 8, 8, 8, 10, 8, 10, 8, 9, 8, 8, 10, 10, 7, 9, 10, 10, 7, 7, 9, 9, 8, 9, 8, 11, 9, 9, 10, 11, 7, 10, 8, 8, 9, 9, 9, 11, 7, 11, 7, 10, 9, 9, 9, 11, 8, 8, 9, 9, 8, 8, 10, 10, 9, 10, 9, 10, 10, 10, 11, 12 }; 
		n &= 63;
		*cc = prec_c[n];
		*d = prec_d[n];
		return prec_a[n];

	}
	if (k == 7) {
		static uint16_t prec_d[128] = {0, 175, 74, 47, 148, 121, 94, 73, 56, 373, 38, 133, 20, 11, 146, 65, 112, 53, 206, 125, 76, 67, 46, 19, 40, 151, 22, 925, 292, 211, 130, 697, 32, 143, 106, 79, 412, 331, 250, 251, 8, 1087, 2, 101, 92, 65, 38, 227, 80, 85, 58, 31, 44, 35, 574, 331, 152, 119, 98, 85, 260, 179, 422, 665, 64, 37, 286, 205, 28, 19, 10, 283, 232, 103, 178, 23, 124, 97, 502, 259, 16, 223, 682, 439, 4, 1, 202, 121, 184, 29, 130, 287, 76, 49, 454, 211, 160, 493, 14, 5, 116, 89, 62, 817, 88, 341, 70, 367, 340, 259, 178, 745, 304, 7, 238, 157, 196, 169, 170, 89, 520, 421, 358, 319, 844, 601, 1330, 2059 };
static uint16_t prec_c[128] = {1, 81, 27, 27, 27, 27, 27, 81, 9, 243, 9, 81, 9, 9, 81, 81, 9, 27, 81, 81, 9, 9, 27, 27, 9, 81, 9, 729, 81, 81, 81, 729, 3, 81, 27, 27, 81, 81, 81, 243, 3, 729, 3, 81, 27, 27, 27, 243, 9, 27, 27, 27, 9, 9, 243, 243, 27, 81, 27, 81, 81, 81, 243, 729, 3, 27, 81, 81, 9, 9, 9, 243, 27, 81, 27, 27, 27, 27, 243, 243, 3, 81, 243, 243, 3, 3, 81, 81, 27, 27, 27, 243, 27, 27, 243, 243, 9, 243, 9, 9, 27, 27, 27, 729, 9, 243, 9, 243, 81, 81, 81, 729, 27, 9, 81, 81, 27, 27, 81, 81, 81, 243, 81, 243, 243, 243, 729, 2187 };
static uint8_t prec_a[128] = {7, 11, 10, 10, 10, 10, 10, 11, 9, 12, 9, 11, 9, 9, 11, 11, 9, 10, 11, 11, 9, 9, 10, 10, 9, 11, 9, 13, 11, 11, 11, 13, 8, 11, 10, 10, 11, 11, 11, 12, 8, 13, 8, 11, 10, 10, 10, 12, 9, 10, 10, 10, 9, 9, 12, 12, 10, 11, 10, 11, 11, 11, 12, 13, 8, 10, 11, 11, 9, 9, 9, 12, 10, 11, 10, 10, 10, 10, 12, 12, 8, 11, 12, 12, 8, 8, 11, 11, 10, 10, 10, 12, 10, 10, 12, 12, 9, 12, 9, 9, 10, 10, 10, 13, 9, 12, 9, 12, 11, 11, 11, 13, 10, 9, 11, 11, 10, 10, 11, 11, 11, 12, 11, 12, 12, 12, 13, 14 };
		n &= 127;
		*cc = prec_c[n];
		*d = prec_d[n];
		return prec_a[n];
	}
#endif
	uint64_t k1 = (k >> 1), k2;
#if STRIP_MSBITS
	cc2 = bzh64(n, k1);
#else
	cc2 = n;
#endif
	ans = helper64(cc2, &d1, &cc1, k1);
	// x->(3^c1*x+d1)/2^k1
	n = (n * cc1 + d1) >> k1;

	k2 = k - k1;
#if STRIP_MSBITS
	n = bzh64(n, k2);
#endif
	ans += helper64(n, &d2, &cc2, k2);
	// (3^c2*(3^c1*x+d1)/2^k1+d2)/2^k2=(3^(c1+c2)*x+(3^c2*d1+d2*2^k1))/2^k
	*d = (d1 * cc2) + (d2 << k1);
	*cc = cc1 * cc2;

	return ans;
}

// n---> (3^c*n+d)/2^k
static uint64_t helper128(uint128_t n, uint64_t * d, uint64_t * cc, uint64_t k)
{
	if (k < THRESHOLD64) {
		// stop using 128 bits numbers and use a 64-bit recursive version
		// when intermediate numbers 
		// are small enough to fit 64 bits registers.

		if (k == 1) {
			// harcoded small values (branch-free)
			uint64_t tn = n & 1;
			*cc = 1 + 2 * tn;
			*d = tn;
			return 1 + tn;
		}

		if (k == 0) {
			*cc = 1;
			*d = 0;
			return 0;
		}

		return helper64((uint64_t) n, d, cc, k);
	}

	uint64_t ans;
	uint64_t d1, d2, cc1, cc2;
	uint128_t r;

	uint64_t k1 = (k >> 1);
#if STRIP_MSBITS
	r = bzh128(n, k1);
#else
	r = n;
#endif
	ans = helper128(r, &d1, &cc1, k1);
	// x->(3^c1*x+d1)/2^k1
	n = (n * cc1 + d1) >> k1;

	uint64_t k2 = k - k1;
#if STRIP_MSBITS
	n = bzh128(n, k2);
#endif
	ans += helper128(n, &d2, &cc2, k2);
	// (3^c2*(3^c1*x+d1)/2^k1+d2)/2^k2=(3^(c1+c2)*x+(3^c2*d1+d2*2^k1))/2^k
	*d = (d1 * cc2) + (d2 << k1);
	*cc = cc1 * cc2;
	/*
	uint64_t hi = *d >> 64;
	if (hi) { printf("ovfl1\n"); exit(1); }
	hi = *cc >> 64;
	if (hi) { printf("ovfl2\n"); exit(1); }
	*/

	return ans;
}

// n---> (3^c*n+d)/2^k
static uint64_t helper(mpz_t n, mpz_t d, mpz_t cc, uint64_t k)
{
	if (k < THRESHOLD128) {
		if (k < THRESHOLD64) {
			// stop using gmp and use a 64-bit recursive version
			// when intermediate numbers 
			// are small enough to fit 64 bits registers.

			if (k == 1) {
				// harcoded small values (branch-free)
				uint64_t tn = mpz_get_ui(n) & 1;
				mpz_set_ui(cc, 1 + 2 * tn);
				mpz_set_ui(d, tn);
				return 1 + tn;
			}

			if (k == 0) {
				mpz_set_ui(cc, 1);
				mpz_set_ui(d, 0);
				return 0;
			}

			uint64_t ans, td, tcc;
			ans = helper64(mpz_get_ui(n), &td, &tcc, k);
			mpz_set_ui(d, td);
			mpz_set_ui(cc, tcc);
			return ans;
		} else {
			// stop using gmp and use a 128-bit recursive version
			// when intermediate numbers 
			// are small enough to fit 128 bits structures optimized by the compiler.
			uint64_t ans;
			uint128_t tn = mpz_get_ui128(n);
			uint64_t tcc, td;
			ans = helper128(tn, &td, &tcc, k);
			mpz_set_ui(d, td);
			mpz_set_ui(cc, tcc);
			return ans;
		}
	}

	uint64_t ans;
	mpz_t d2, cc2;
	mpz_init(cc2);
	mpz_init(d2);

	uint64_t k1 = (k >> 1);
	uint64_t k2 = k - k1;

	mpz_fdiv_r_2exp(cc2, n, k1);
	ans = helper(cc2, d, cc, k1);
	// x->(3^c1*x+d1)/2^k1
#if 0
	if (mpz_sizeinbase(cc, 2) + mpz_sizeinbase(n, 2) + 1 > k)
	 {
        		printf("cc %d n %d d %d k %d\n", mpz_sizeinbase(cc, 2) , mpz_sizeinbase(n, 2) , mpz_sizeinbase(d, 2), k);
	 }

#endif
#if 0
	mpz_mul(n, n, cc);
	mpz_add(n, n, d);
	mpz_fdiv_q_2exp(n, n, k1);

	mpz_fdiv_r_2exp(n, n, k2);
#else
	mpz_mullo(n, k, n, cc);
	mpz_add(n, n, d);
	mpz_fdiv_q_2exp(n, n, k1);
	mpz_fdiv_r_2exp(n, n, k2);
#endif
	ans += helper(n, d2, cc2, k2);
	// (3^c2*(3^c1*x+d1)/2^k1+d2)/2^k2=(3^(c1+c2)*x+(3^c2*d1+d2*2^k1))/2^k
	mpz_mul(d, d, cc2);
	mpz_mul_2exp(d2, d2, k1);
	mpz_add(d, d, d2);
	mpz_mul(cc, cc, cc2);

	mpz_clear(cc2);
	mpz_clear(d2);
	return ans;
}

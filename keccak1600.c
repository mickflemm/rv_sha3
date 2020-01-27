// SPDX-License-Identifier: BSD-2-Clause
/*
 * Keccak-f[1600] C Implementation - State permutation
 * Copyright (C) 2019 Nick Kossifidis <mick@ics.forth.gr>
 */

#include "keccak1600.h"

/************************\
* PRE-COMPUTED CONSTANTS *
\************************/

#if 0
/*
 * Pre-computed lane indices for the pi step, where
 * (x,y) = ((0,1) * (2,3)) * (x,y) =>
 * (x,y) = (0 * x + 1 * y, 2 * x + 3 * y)
 */
static const uint8_t pi_lane_idxes[KECCAK_NUM_LANES] =
	{  0, 10, 20,  5, 15,	/* y = 0 */
	  16,  1, 11, 21,  6,	/* y = 1 */
	   7, 17,  2, 12, 22,	/* y = 2 */
	  23,  8, 18,  3, 13,	/* y = 3 */
	  14, 24,  9, 19,  4};	/* y = 4 */
/*x:	   0, 1,  2,  3,  4 */

/*
 * Pre-computed rotation offsets for the rho step,
 * from Keccak Reference, Table 2.1, the final table
 * (after modulo <lane size>) here can be found on
 * https://keccak.team/keccak_specs_summary.html
 */
static const uint8_t rho_offsets[KECCAK_NUM_LANES] =
	{  0,  1, 62, 28, 27,	/* y = 0 */
	  36, 44,  6, 55, 20,	/* y = 1 */
	   3, 10, 43, 25, 39,	/* y = 2 */
	  41, 45, 15, 21,  8,	/* y = 3 */
	  18,  2, 61, 56, 14};	/* y = 4 */
/*x:	   0, 1,  2,  3,  4 */
#endif

/*
 * If we follow the above approach when applying rho and pi, we'll end up
 * needing a copy of the state so that we don't use values that we already
 * modified, wasting up memory for no reason. Instead of working in the
 * order (0, 0) -> (1, 0) -> (2, 0)... shown on the tables above, we can
 * start from point (1, 0) and proceed by following the pi mapping.
 *
 * We'll go from (x,y) to (y, 2*x + 3y) on each step, so it's going to be
 * (1, 0) -> (0, 2) -> (2, 1) -> (1, 2)... This way we won't hit the same
 * lane twice and we can get away by just saving the first row (or as shown
 * in rho_pi below, just the first lane (1,0)). You may also check out
 * keccakf1600_round_unrolled below for easier inspection. The idea comes
 * from section 2.5 of "Keccak implementation overview"
 */
#if !(defined(USE_UNROLLED) || defined(USE_UNROLLED_CF))
static const uint8_t pi_lane_idxes[KECCAK_NUM_LANES - 1] =
	{ 10,  7, 11, 17, 18,
	   3,  5, 16,  8, 21,
	  24,  4, 15, 23, 19,
	  13, 12,  2, 20, 14,
	  22,  9,  6,  1};

static const uint8_t rho_offsets[KECCAK_NUM_LANES - 1] =
	{ 1,  3,  6, 10, 15,
	 21, 28, 36, 45, 55,
	  2, 14, 27, 41, 56,
	  8, 25, 43, 62, 18,
	 39, 61, 20, 44};
#endif

/*
 * Pre-computed round constants for the Iota step, Keccak Reference,
 * Section 1.2. If this takes too much memory we can also use the
 * in-place calculation (see iota below).
 */
#if !defined(__OPTIMIZE_SIZE__) || (defined(USE_UNROLLED) || defined(USE_UNROLLED_CF))
static const lane_t round_constants[KECCAK1600_NUM_ROUNDS] =
    { 0x0000000000000001ULL, 0x0000000000008082ULL, 0x800000000000808aULL,
	0x8000000080008000ULL, 0x000000000000808bULL, 0x0000000080000001ULL,
	0x8000000080008081ULL, 0x8000000000008009ULL, 0x000000000000008aULL,
	0x0000000000000088ULL, 0x0000000080008009ULL, 0x000000008000000aULL,
	0x000000008000808bULL, 0x800000000000008bULL, 0x8000000000008089ULL,
	0x8000000000008003ULL, 0x8000000000008002ULL, 0x8000000000000080ULL,
	0x000000000000800aULL, 0x800000008000000aULL, 0x8000000080008081ULL,
	0x8000000000008080ULL, 0x0000000080000001ULL, 0x8000000080008008ULL };
#endif

/*********\
* HELPERS *
\*********/

/*
 * Lookup table for modulo 5, to speed things up.
 * We know that the largest value we'll hit will
 * be x + 4 where x = 4, so we can get away with
 * 9 bytes.
 */
static const uint8_t mod5[9] = { 0, 1, 2, 3, 4, 0, 1, 2, 3 };

/* Left-rotate a lane */
static inline lane_t rotl_lane(lane_t val, int times)
{
	return (((val) << (times)) |
		((val) >> (KECCAK1600_LANE_BITS - (times))));
}

#if !defined(USE_UNROLLED) && !defined(USE_UNROLLED_CF)

/**********************\
* KECCAK STEP MAPPINGS *
\**********************/

/*
 * Theta step, Keccak Reference, Section 2.3.2
 */
static inline void theta(lane_t * A)
{
	unsigned int x = 0;
	unsigned int y_offset = 0;
	lane_t C[5] = { 0 };
	lane_t D = 0;

	/* Compute the parity of the columns */
	for (x = 0; x < KECCAL_NUM_COLS; x++)
		/* (x, 0) ^ (x, 1) ^ (x, 2) ^ (x, 3) ^ (x, 4) */
		C[x] = A[x] ^ A[x + 5] ^ A[x + 10] ^ A[x + 15] ^ A[x + 20];

	for (x = 0; x < KECCAL_NUM_COLS; x++) {
		/* Compute theta for this column */
		D = C[mod5[x + 4]] ^ rotl_lane(C[mod5[x + 1]], 1);

		/* Apply theta to each row of this slice (x-y plane) */
		for (y_offset = 0; y_offset < KECCAK_NUM_LANES; y_offset += 5)
			A[x + y_offset] ^= D;
	}
}

/*
 * Rho and pi steps, Keccak Reference, Sections 2.3.3 and 2.3.4
 * These two steps are easily combined, rho step is the rotation
 * of A[x,y] by rho_offsets[x, y] and pi step is the mapping
 * of A[x,y] -> A[y, 2x + 3y]. So we do:
 * A[y, 2x + 3y] = rotl64(A[x,y], rho_offsets[x,y])
 */
static inline void rho_pi(lane_t * A)
{
	unsigned int i = 0;
	lane_t start = 0;
	lane_t previous = 0;

	/* If we move forward we'll need to save the first row
	 * as mentioned above, however if we move backwards we
	 * only need to store the value of (1,0) to be used when
	 * we reach the begining. */
	start = A[1];		/* (1, 0) */
	for (i = KECCAK_NUM_LANES - 2; i > 0; i--) {
		previous = A[pi_lane_idxes[i - 1]];
		A[pi_lane_idxes[i]] = rotl_lane(previous, rho_offsets[i]);
	}
	/* Reached (0,2), move to (1,0) */
	A[pi_lane_idxes[0]] = rotl_lane(start, rho_offsets[0]);
}

/*
 * Chi step, Keccak Reference, Section 2.3.1
 */
static inline void chi(lane_t * A)
{
	unsigned int x = 0;
	unsigned int y_offset = 0;
	lane_t a[5] = { 0 };

	for (y_offset = 0; y_offset < KECCAK_NUM_LANES; y_offset += 5) {
		/* Make a copy of the plane and then
		 * apply chi on it */
		for (x = 0; x < KECCAL_NUM_COLS; x++)
			a[x] = A[x + y_offset];

		for (x = 0; x < KECCAL_NUM_COLS; x++)
			A[x + y_offset] ^= (~a[mod5[x + 1]] & a[mod5[x + 2]]);
	}

}

/*
 * Iota step, Keccak Reference, Section 2.3.5
 */
#ifndef __OPTIMIZE_SIZE__
static inline void iota(lane_t * A, int round)
{
	A[0] ^= round_constants[round];
}
#else
static __attribute__ ((noinline))
void iota(lane_t * A, uint8_t * reg)
{
	int i = 0;
	unsigned int bit_i = 0;

	for (i = 0; i < 7; i++) {
		/* bit position 2^i - 1 */
		bit_i = (1 << i) - 1;
		if ((*reg) & 0x01)
			A[0] ^= (lane_t) 1 << bit_i;
		if ((*reg) & 0x80)
			/* LFSR Register using
			 * Primitive polynomial over
			 * GF(2): x^8+x^6+x^5+x^4+1
			 */
			(*reg) = ((*reg) << 1) ^ 0x71;
		else
			(*reg) <<= 1;
	}
}
#endif


/**********************************************\
* KECCAK-F1600 STATE PERMUTATION / ENTRY POINT *
\**********************************************/

void
keccakf1600_state_permute(k1600_state_t * st)
{
#ifdef __OPTIMIZE_SIZE__
	uint8_t lfsr_reg = 0x1;
#endif
	int i = 0;
	for (i = 0; i < KECCAK1600_NUM_ROUNDS; i++) {
		theta(st->A);
		rho_pi(st->A);
		chi(st->A);
#ifdef __OPTIMIZE_SIZE__
		iota(st->A, &lfsr_reg);
#else
		iota(st->A, i);
#endif
	}

}

#elif defined(USE_UNROLLED)

/*
 * This is an unrolled version of the above
 * with some further cleanups to see how well
 * it competes with the compiler's optimizations.
 * It also provides a better view of what's going on.
 */
static void
keccakf1600_round_unrolled(lane_t * A, int r_idx)
{
	register lane_t C[5] = { 0 };
	register lane_t T = 0;

	/* Compute the parity of the columns */

	C[0] = A[0] ^ A[5] ^ A[10] ^ A[15] ^ A[20];
	C[1] = A[1] ^ A[6] ^ A[11] ^ A[16] ^ A[21];
	C[2] = A[2] ^ A[7] ^ A[12] ^ A[17] ^ A[22];
	C[3] = A[3] ^ A[8] ^ A[13] ^ A[18] ^ A[23];
	C[4] = A[4] ^ A[9] ^ A[14] ^ A[19] ^ A[24];

	/* Compute theta for each column */

	T = C[4];
	C[4] ^= rotl_lane(C[1], 1);	// D[0]
	C[1] ^= rotl_lane(C[3], 1);	// D[2]
	C[3] ^= rotl_lane(C[0], 1);	// D[4]
	C[0] ^= rotl_lane(C[2], 1);	// D[1]
	C[2] ^= rotl_lane(T, 1);	// D[3]

	/* Apply theta-rho-pi in-place */

	A[0] ^= C[4];

	T = A[1];
	A[1]  = rotl_lane(A[6]  ^ C[0], 44);
	A[6]  = rotl_lane(A[9]  ^ C[3], 20);
	A[9]  = rotl_lane(A[22] ^ C[1], 61);
	A[22] = rotl_lane(A[14] ^ C[3], 39);
	A[14] = rotl_lane(A[20] ^ C[4], 18);
	A[20] = rotl_lane(A[2]  ^ C[1], 62);
	A[2]  = rotl_lane(A[12] ^ C[1], 43);
	A[12] = rotl_lane(A[13] ^ C[2], 25);
	A[13] = rotl_lane(A[19] ^ C[3], 8);
	A[19] = rotl_lane(A[23] ^ C[2], 56);
	A[23] = rotl_lane(A[15] ^ C[4], 41);
	A[15] = rotl_lane(A[4]  ^ C[3], 27);
	A[4]  = rotl_lane(A[24] ^ C[3], 14);
	A[24] = rotl_lane(A[21] ^ C[0], 2);
	A[21] = rotl_lane(A[8]  ^ C[2], 55);
	A[8]  = rotl_lane(A[16] ^ C[0], 45);
	A[16] = rotl_lane(A[5]  ^ C[4], 36);
	A[5]  = rotl_lane(A[3]  ^ C[2], 28);
	A[3]  = rotl_lane(A[18] ^ C[2], 21);
	A[18] = rotl_lane(A[17] ^ C[1], 15);
	A[17] = rotl_lane(A[11] ^ C[0], 10);
	A[11] = rotl_lane(A[7]  ^ C[1], 6);
	A[7]  = rotl_lane(A[10] ^ C[4], 3);
	A[10] = rotl_lane(T ^ C[0], 1);

	/* Apply chi on each plane */

	C[0] = A[0];
	C[1] = A[1];

	A[0] ^= (~A[1] & A[2]);
	A[1] ^= (~A[2] & A[3]);
	A[2] ^= (~A[3] & A[4]);
	A[3] ^= (~A[4] & C[0]);
	A[4] ^= (~C[0] & C[1]);

	C[0] = A[5];
	C[1] = A[6];

	A[5] ^= (~A[6] & A[7]);
	A[6] ^= (~A[7] & A[8]);
	A[7] ^= (~A[8] & A[9]);
	A[8] ^= (~A[9] & C[0]);
	A[9] ^= (~C[0] & C[1]);

	C[0] = A[10];
	C[1] = A[11];

	A[10] ^= (~A[11] & A[12]);
	A[11] ^= (~A[12] & A[13]);
	A[12] ^= (~A[13] & A[14]);
	A[13] ^= (~A[14] & C[0]);
	A[14] ^= (~C[0] & C[1]);

	C[0] = A[15];
	C[1] = A[16];

	A[15] ^= (~A[16] & A[17]);
	A[16] ^= (~A[17] & A[18]);
	A[17] ^= (~A[18] & A[19]);
	A[18] ^= (~A[19] & C[0]);
	A[19] ^= (~C[0] & C[1]);

	C[0] = A[20];
	C[1] = A[21];

	A[20] ^= (~A[21] & A[22]);
	A[21] ^= (~A[22] & A[23]);
	A[22] ^= (~A[23] & A[24]);
	A[23] ^= (~A[24] & C[0]);
	A[24] ^= (~C[0] & C[1]);

	/* Apply iota */

	A[0] ^= round_constants[r_idx];
}

void
keccakf1600_state_permute(k1600_state_t * st)
{
	int i = 0;
	for (i = 0; i < KECCAK1600_NUM_ROUNDS; i++)
		keccakf1600_round_unrolled(st->A, i);
}

#elif defined(USE_UNROLLED_CF)

/*
 * On the above approach when modifying A in-place, we
 * go back and forth accessing lanes that are often
 * further away than a cache line (e.g. 24 -> 4), resulting
 * cache line misses. Here instead of modifying A in-place,
 * we switch between A and an intermediate state so that
 * we can load from one and store to the other. Since
 * the number of rounds is even we will end up with A
 * at the end. We also access lanes that are no further
 * than 6 lanes away each time. This way we get fewer
 * cache-line misses and also allow OoO CPUs to be more
 * efficient.
 *
 * This is the default approach in OpenSSL (KECCAK_2X)
 */
static void
keccakf1600_round_unrolled_cf(lane_t * A, lane_t * N, int r_idx)
{
	register lane_t C[5] = { 0 };
	register lane_t T[5] = { 0 };

	/* Compute the parity of the columns */

	C[0] = A[0] ^ A[5] ^ A[10] ^ A[15] ^ A[20];
	C[1] = A[1] ^ A[6] ^ A[11] ^ A[16] ^ A[21];
	C[2] = A[2] ^ A[7] ^ A[12] ^ A[17] ^ A[22];
	C[3] = A[3] ^ A[8] ^ A[13] ^ A[18] ^ A[23];
	C[4] = A[4] ^ A[9] ^ A[14] ^ A[19] ^ A[24];

	/* Compute theta for each column */

	T[0] = C[4];
	C[4] ^= rotl_lane(C[1], 1);	// D[0]
	C[1] ^= rotl_lane(C[3], 1);	// D[2]
	C[3] ^= rotl_lane(C[0], 1);	// D[4]
	C[0] ^= rotl_lane(C[2], 1);	// D[1]
	C[2] ^= rotl_lane(T[0], 1);	// D[3]

	/* Apply theta-rho-pi */

	T[0] = A[0] ^ C[4];
	T[1] = rotl_lane(A[6]  ^ C[0], 44);
	T[2] = rotl_lane(A[12] ^ C[1], 43);
	T[3] = rotl_lane(A[18] ^ C[2], 21);
	T[4] = rotl_lane(A[24] ^ C[3], 14);

	/* Apply chi on the 1st plane*/

	/* Also apply iota since we are here */
	N[0] = T[0] ^ (~T[1] & T[2]) ^ round_constants[r_idx];
	N[1] = T[1] ^ (~T[2] & T[3]);
	N[2] = T[2] ^ (~T[3] & T[4]);
	N[3] = T[3] ^ (~T[4] & T[0]);
	N[4] = T[4] ^ (~T[0] & T[1]);


	T[0] = rotl_lane(A[3]  ^ C[2], 28);
	T[1] = rotl_lane(A[9]  ^ C[3], 20);
	T[2] = rotl_lane(A[10] ^ C[4], 3);
	T[3] = rotl_lane(A[16] ^ C[0], 45);
	T[4] = rotl_lane(A[22] ^ C[1], 61);

	N[5] = T[0] ^ (~T[1] & T[2]);
	N[6] = T[1] ^ (~T[2] & T[3]);
	N[7] = T[2] ^ (~T[3] & T[4]);
	N[8] = T[3] ^ (~T[4] & T[0]);
	N[9] = T[4] ^ (~T[0] & T[1]);


	T[0] = rotl_lane(A[1]  ^ C[0], 1);
	T[1] = rotl_lane(A[7]  ^ C[1], 6);
	T[2] = rotl_lane(A[13] ^ C[2], 25);
	T[3] = rotl_lane(A[19] ^ C[3], 8);
	T[4] = rotl_lane(A[20] ^ C[4], 18);

	N[10] = T[0] ^ (~T[1] & T[2]);
	N[11] = T[1] ^ (~T[2] & T[3]);
	N[12] = T[2] ^ (~T[3] & T[4]);
	N[13] = T[3] ^ (~T[4] & T[0]);
	N[14] = T[4] ^ (~T[0] & T[1]);


	T[0] = rotl_lane(A[4]  ^ C[3], 27);
	T[1] = rotl_lane(A[5]  ^ C[4], 36);
	T[2] = rotl_lane(A[11] ^ C[0], 10);
	T[3] = rotl_lane(A[17] ^ C[1], 15);
	T[4] = rotl_lane(A[23] ^ C[2], 56);

	N[15] = T[0] ^ (~T[1] & T[2]);
	N[16] = T[1] ^ (~T[2] & T[3]);
	N[17] = T[2] ^ (~T[3] & T[4]);
	N[18] = T[3] ^ (~T[4] & T[0]);
	N[19] = T[4] ^ (~T[0] & T[1]);


	T[0] = rotl_lane(A[2]  ^ C[1], 62);
	T[1] = rotl_lane(A[8]  ^ C[2], 55);
	T[2] = rotl_lane(A[14] ^ C[3], 39);
	T[3] = rotl_lane(A[15] ^ C[4], 41);
	T[4] = rotl_lane(A[21] ^ C[0], 2);

	N[20] = T[0] ^ (~T[1] & T[2]);
	N[21] = T[1] ^ (~T[2] & T[3]);
	N[22] = T[2] ^ (~T[3] & T[4]);
	N[23] = T[3] ^ (~T[4] & T[0]);
	N[24] = T[4] ^ (~T[0] & T[1]);
}

void
keccakf1600_state_permute(k1600_state_t * st)
{
	lane_t T[KECCAK_NUM_LANES] = { 0 };
	int i = 0;

	for (i = 0; i < KECCAK1600_NUM_ROUNDS; i += 2) {
		keccakf1600_round_unrolled_cf(st->A, T, i);
		keccakf1600_round_unrolled_cf(T, st->A, i + 1);
	}
}

#endif

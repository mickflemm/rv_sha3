// SPDX-License-Identifier: BSD-2-Clause
/*
 * Keccak-f[1600] C Implementation - State permutation
 * Copyright (C) 2019 Nick Kossifidis <mick@ics.forth.gr>
 */

#include "keccak1600.h"


/**********************\
* KECCAK STEP MAPPINGS *
\**********************/

/*
 * Theta step, Keccak Reference, Section 2.3.2
 *
 * This is a linear transformation that acts as the main diffusion layer.
 * It gets the sum of parity bits of columns (x - 1, *, z),  (x + 1, *, z - 1)
 * and adds them to every row (for each y) of the slice (x-y plane).
 */
static inline void theta(lane_t *A)
{
	uint_fast8_t x = 0;
	uint_fast8_t y_offset = 0;
	lane_t C[KECCAK_NUM_COLS] = { 0 };
	lane_t D[KECCAK_NUM_COLS] = { 0 };

	/* Compute the sum of parity bits of each column */
	for (x = 0; x < KECCAK_NUM_COLS; x++)
		/* (x, 0) ^ (x, 1) ^ (x, 2) ^ (x, 3) ^ (x, 4) */
		C[x] = A[x] ^ A[x + 5] ^ A[x + 10] ^ A[x + 15] ^ A[x + 20];

	/* Compute D for each slice:
	 * C of column ((x - 1) % 5, *, z) ^ C of column ((x + 1) % 5, *, z -1)
	 * (we unroll this so no no need for % 5, rotl_lane by 1 is for the z-1 part)
	 */
	D[0] = C[4] ^ rotl_lane(C[1], 1);
	D[1] = C[0] ^ rotl_lane(C[2], 1);
	D[2] = C[1] ^ rotl_lane(C[3], 1);
	D[3] = C[2] ^ rotl_lane(C[4], 1);
	D[4] = C[3] ^ rotl_lane(C[0], 1);

	for (x = 0; x < KECCAK_NUM_COLS; x++) {
		/* Apply D to each row of this slice */
		for (y_offset = 0; y_offset < KECCAK_NUM_LANES; y_offset += 5)
			A[x + y_offset] ^= D[x];
	}
}

/*
 * Pi step, Section 2.3.3
 *
 * This step is a linear transformation of (x,y) using the
 * matrix: |0,1|, which results (x,y) -> (y, 2*x + 3*y)
 *         |2,3|
 * it's aimed at intra-slice dispersion / long-term diffusion.
 *
 * To avoid calculating (y, 2x + 3y) for each (x,y), we can pre-compute
 * the mappings using A indices, so for example (1,0) which is A[1]
 * will become (0, 2*1 + 3*0) = (0,2) = A[0 + 5*2] = A[10], so we can
 * create an array to store the mapped indices like this:
 * Y
 * 0 |00, 10, 20, 05, 15,
 * 1 |16, 01, 11, 21, 06,
 * 2 |07, 17, 02, 12, 22,
 * 3 |23, 08, 18, 03, 13,
 * 4 |14, 24, 09, 19, 04
 *   -------------------
 * X   0   1   2   3   4
 * (note the diagonals, that's a visual representation of the pi mapping)
 * 
 * With this approach we 'll work in a sequence from (0,0) to (4,4), and
 * end up using lanes we 've already modified, for example after A[2] becomes
 * A[20], we 'won't be able to set A[12] to A[2] unless we keep a copy of
 * the whole state. A better approach is to take advantage of the pi mapping
 * and work in that order, from (x,y) to (y, 2x+3y) on each step, so instead of
 * working in the same order as A, from (1,0) -> (2,0) -> (3,0) etc, it's going
 * to be (1, 0) -> (0, 2) -> (2, 1) -> (1, 2)... This way we won't hit the same
 * lane twice, and only depend on lanes we haven't modified and the first row
 * (x,0) that we 'll reach again by the end of the sequence. An even better
 * approach is while following the pi mapping, to start from the end instead
 * of (1,0) and work backwards, so that we avoid modifying the first row, and
 * only store the first lane for the final step (0,2) -> (1,0). The idea comes
 * from section 2.5 of "Keccak implementation overview".
 *
 * With the above in mind, here is the final array of indices for the Pi mapping
 * excluding 0 (since pi mapping doesn't modify 0,0).
 */
static const uint_fast8_t pi_lane_idxes[KECCAK_NUM_LANES - 1] =
	{ 10,  7, 11, 17, 18,
	   3,  5, 16,  8, 21,
	  24,  4, 15, 23, 19,
	  13, 12,  2, 20, 14,
	  22,  9,  6,  1};

/*
 * Rho step, Section 2.3.4
 *
 * This is another linear transformation for improving inter-slice
 * dispersion. It's rotating each lane by a constant coming from the
 * sequence i*(i+1)/2 = (0, 1, 3, 6, 10, 15, 21, 28, 36, 45, 55, 66...),
 * modulo the lane size (64). This sequence guarantees that all constants
 * are different and not present more than once. The rotation is applied to
 * (x,y), multiplied by the same transformation matrix as pi, here are the
 * constants for each (x,y) on A:
 *
 * Y
 * 0 |00, 01, 62, 28, 27,
 * 1 |36, 44, 06, 55, 20,
 * 2 |03, 10, 43, 25, 39,
 * 3 |41, 45, 15, 21, 08,
 * 4 |18, 02, 61, 56, 14
 *   -------------------
 * X   0   1   2   3   4
 * 
 * (comes from https://keccak.team/keccak_specs_summary.html)
 *
 * However since we follow the pi mapping anyway, we can apply rho together with pi,
 * and this array just becomes the same as the sequence, we keep it here for convenience:
 */
 #if !defined(__OPTIMIZE_SIZE__)
static const uint_fast8_t rho_offsets[KECCAK_NUM_LANES - 1] =
	{ 1,  3,  6, 10, 15,
	 21, 28, 36, 45, 55,
	  2, 14, 27, 41, 56,
	  8, 25, 43, 62, 18,
	 39, 61, 20, 44};
static inline int get_rho_for_idx(int idx)
{
	return rho_offsets[idx];
}
#else
 /* This will return the rotation constants for each index
 * of pi_lane_idxes, since it's the sequence (i*(i+1) / 2)
 * mod 64 (when following the pi mapping), it's simple/fast
 * to calculate it and is still smaller than the 23byte lookup table. */
static inline int get_rho_for_idx(int idx)
{
	int i = idx + 1; /* Start from 1 since 0,0 is ignored */
	return ((i * (i + 1)) >> 1) & 0x3F;
}
#endif

static inline void rho_pi(lane_t *A)
{
	uint_fast8_t i = 0;
	lane_t first = A[1]; /* Save (1,0) for the last step */
	for (i = KECCAK_NUM_LANES - 2; i > 0; i--) {
		lane_t next = A[pi_lane_idxes[i - 1]];
		A[pi_lane_idxes[i]] = rotl_lane(next, get_rho_for_idx(i));
	}
	/* Reached (0,2), move to (1,0) */
	A[pi_lane_idxes[0]] = rotl_lane(first, get_rho_for_idx(0));
}

/*
 * Chi step, Keccak Reference, Section 2.3.1
 *
 * This is the only non-linear transformation, without it a Keccak
 * round would be a linear function. The formula is:
 * A[(x,y)] = A[(x,y)] ^ (~A[((x + 1) % 5, y)] & A[((x + 2) %5, y)])
 *
 * Note: There is a trick here to avoid most NOTs, but it looks better
 * in the fully-unrolled version.
 */
static inline void chi(lane_t *A)
{
	uint_fast8_t y_offset = 0;
	lane_t a[KECCAK_NUM_COLS] = { 0 };

	for (y_offset = 0; y_offset < KECCAK_NUM_LANES; y_offset += 5) {
		uint_fast8_t x = 0;
		/* Make a copy of the plane and then
		 * apply chi on it */
		for (x = 0; x < KECCAK_NUM_COLS; x++)
			a[x] = A[x + y_offset];

		A[0 + y_offset] ^= (~a[1] & a[2]);
		A[1 + y_offset] ^= (~a[2] & a[3]);
		A[2 + y_offset] ^= (~a[3] & a[4]);
		A[3 + y_offset] ^= (~a[4] & a[0]);
		A[4 + y_offset] ^= (~a[0] & a[1]);
	}
}

/*
 * Iota step, Keccak Reference, Section 2.3.5
 *
 * This step is another linear transformation, aimed at breaking
 * symmetry across rounds. It defines a set of constants to be xored
 * to the first lane on each round. The formula for calculating
 * the constants, based on an LFSR is defined on section 1.2, and
 * a thing to notice is that there are only 7bits that are being
 * modified each time, all the rest are zeroes, so the round constants
 * can be easily compressed. We'll use that property for the
 * size-optimized implementation. Here we'll just store the round
 * constants to be xored as the authors intended, the table comes from
 * https://keccak.team/keccak_specs_summary.html
 */
 #if !defined(__OPTIMIZE_SIZE__)
 static const lane_t round_constants[KECCAK1600_NUM_ROUNDS] =
 	{ 0x0000000000000001ULL, 0x0000000000008082ULL, 0x800000000000808aULL,
	  0x8000000080008000ULL, 0x000000000000808bULL, 0x0000000080000001ULL,
	  0x8000000080008081ULL, 0x8000000000008009ULL, 0x000000000000008aULL,
	  0x0000000000000088ULL, 0x0000000080008009ULL, 0x000000008000000aULL,
	  0x000000008000808bULL, 0x800000000000008bULL, 0x8000000000008089ULL,
	  0x8000000000008003ULL, 0x8000000000008002ULL, 0x8000000000000080ULL,
	  0x000000000000800aULL, 0x800000008000000aULL, 0x8000000080008081ULL,
	  0x8000000000008080ULL, 0x0000000080000001ULL, 0x8000000080008008ULL };

static inline void iota(lane_t *A, unsigned int round)
{
	A[0] ^= round_constants[round];
}
#else
/* Storing round constants in 64bit format takes a lot of space, calculating them
 * using the LFSR on each round is on the other hand quite slow. Here we take advantage
 * of the fact that rc constants only use bits 2^i - 1 (so 0, 1, 3, 7, 15, 31, 63), so
 * we can map those bits to a 7bit set, and use a byte instead of 8bytes for each rc. */
static const uint8_t rc_compressed[KECCAK1600_NUM_ROUNDS] =
	{ 0x01, 0x1A, 0x5E, 0x70, 0x1F, 0x21, 0x79, 0x55,
	  0x0E, 0x0C, 0x35, 0x26, 0x3F, 0x4F, 0x5D, 0x53,
	  0x52, 0x48, 0x16, 0x66, 0x79, 0x58, 0x21, 0x74 };

static inline void iota(lane_t *A, unsigned int round)
{
	uint64_t rc = 0;
	uint8_t rc_comp = rc_compressed[round];
	int i = 0;
	for (i = 0; i < 7; i++) {
		/* bit position on rc (2^i - 1) */
		unsigned int bit_i = (1 << i) - 1;
		/* bit position on compressed rc */
		unsigned int bit_c = 1 << i;
		if(rc_comp & bit_c)
                        rc |= 1ULL << bit_i;
        }
	A[0] ^= rc;
}
#endif

/**********************************************\
* KECCAK-F1600 STATE PERMUTATION / ENTRY POINT *
\**********************************************/

void
keccakf1600_state_permute_ref(k1600_state_t *st)
{
	int i = 0;
	for (i = 0; i < KECCAK1600_NUM_ROUNDS; i++) {
		theta(st->A);
		rho_pi(st->A);
		chi(st->A);
		iota(st->A, i);
	}
}
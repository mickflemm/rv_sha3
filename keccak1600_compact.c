// SPDX-License-Identifier: BSD-2-Clause
/*
 * Keccak-f[1600] C Implementation - State permutation
 * Copyright (C) 2024 Nick Kossifidis <mick@ics.forth.gr>
 */

#include "keccak1600.h"

/*
 * This is the original implementation, focused on
 * memory / storage constrained environments.
 */

static inline void theta(lane_t *A)
{
	uint_fast8_t x = 0;
	uint_fast8_t y_offset = 0;
	lane_t C[KECCAK_NUM_COLS] = { 0 };

	/* Compute the sum of parity bits of each column */
	for (x = 0; x < KECCAK_NUM_COLS; x++)
		/* (x, 0) ^ (x, 1) ^ (x, 2) ^ (x, 3) ^ (x, 4) */
		C[x] = A[x] ^ A[x + 5] ^ A[x + 10] ^ A[x + 15] ^ A[x + 20];

	for (x = 0; x < KECCAK_NUM_COLS; x++) {
		/* Compute D for this column */
		lane_t D = C[(x + 4) % 5] ^ rotl_lane(C[(x + 1) % 5], 1);

		/* Apply D to each row of this slice */
		for (y_offset = 0; y_offset < KECCAK_NUM_LANES; y_offset += 5)
			A[x + y_offset] ^= D;
	}
}

/* Converting from index to x,y coords, and calculating 
 * 2x + 3y etc takes more space than 24bytes and it's also
 * much slower, even worse if we want to do the in-place
 * approach, work backwards etc. This should fit in a cache
 * line. */
static const uint8_t pi_lane_idxes[KECCAK_NUM_LANES - 1] =
	{ 10,  7, 11, 17, 18,
	   3,  5, 16,  8, 21,
	  24,  4, 15, 23, 19,
	  13, 12,  2, 20, 14,
	  22,  9,  6,  1};

/* This will return the rotation constants for each index
 * of pi_lane_idxes, since it's the sequence (i*(i+1) / 2)
 * mod 64 (when following the pi mapping), it's simple/fast
 * to calculate it and takes in worst case 5 + 32bit instructions
 * (20 bytes) which is still smaller than the 23byte lookup table. */
static inline int get_rho_for_idx(int idx)
{
	int i = idx + 1; /* Start from 1 since 0,0 is ignored */
	return ((i * (i + 1)) >> 1) & 0x3F;
}

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

static inline void chi(lane_t * A)
{
	unsigned int x = 0;
	unsigned int y_offset = 0;
	lane_t a[5] = { 0 };

	for (y_offset = 0; y_offset < KECCAK_NUM_LANES; y_offset += 5) {
		/* Make a copy of the plane and then
		 * apply chi on it */
		for (x = 0; x < KECCAK_NUM_COLS; x++)
			a[x] = A[x + y_offset];

		for (x = 0; x < KECCAK_NUM_COLS; x++)
			A[x + y_offset] ^= (~a[(x + 1) % 5] & a[(x + 2) % 5]);
	}
}

/* Storing round constants in 64bit format takes a lot of space, calculating them
 * using the LFSR on each round is on the other hand quite slow. Here we take advantage
 * of the fact that rc constants only use bits 2^i - 1 (so 0, 1, 3, 7, 15, 31, 63), so
 * we can map those bits to a 7bit set, and use a byte instead of 8bytes for each rc.
 * In case we use 32bit instructions the 24 bytes this array takes is 3 instructions and
 * if we used the LFSR approach we'd need more, so this approach takes less space. */
uint8_t rc_compressed[KECCAK1600_NUM_ROUNDS] =
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

void
keccakf1600_state_permute_compact(k1600_state_t *st)
{
	int i = 0;
	for (i = 0; i < KECCAK1600_NUM_ROUNDS; i++) {
		theta(st->A);
		rho_pi(st->A);
		chi(st->A);
		iota(st->A, i);
	}
}
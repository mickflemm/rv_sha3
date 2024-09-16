// SPDX-License-Identifier: BSD-2-Clause
/*
 * Keccak-f[1600] C Implementation - State permutation
 * Copyright (C) 2019 Nick Kossifidis <mick@ics.forth.gr>
 */

#include "keccak1600.h"

/*
 * When working in-place, we go back and forth accessing lanes
 * that are often further away than a cache line (e.g. 24 -> 4),
 * resulting cache line misses. 
 *
 * Here instead of modifying A in-place, we switch between
 * A and an intermediate state N so that we can load from one
 * and store to the other. Since the number of rounds is even
 * we will end up with A at the end. We also access lanes that
 * are no further than 6 lanes away each time. This way we get
 * fewer cache-line misses and also allow OoO CPUs to be more
 * efficient.
 */

static const lane_t round_constants[KECCAK1600_NUM_ROUNDS] =
	{ 0x0000000000000001ULL, 0x0000000000008082ULL, 0x800000000000808aULL,
	  0x8000000080008000ULL, 0x000000000000808bULL, 0x0000000080000001ULL,
	  0x8000000080008081ULL, 0x8000000000008009ULL, 0x000000000000008aULL,
	  0x0000000000000088ULL, 0x0000000080008009ULL, 0x000000008000000aULL,
	  0x000000008000808bULL, 0x800000000000008bULL, 0x8000000000008089ULL,
	  0x8000000000008003ULL, 0x8000000000008002ULL, 0x8000000000000080ULL,
	  0x000000000000800aULL, 0x800000008000000aULL, 0x8000000080008081ULL,
	  0x8000000000008080ULL, 0x0000000080000001ULL, 0x8000000080008008ULL };

static inline void
keccakf1600_round_intermediate_unrolled(lane_t *A, lane_t *N, int r_idx)
{
	lane_t C[5] = { 0 };
	lane_t D[5] = { 0 };
	lane_t T[5] = { 0 };

	/* Compute parity of columns */
	C[0] = A[0] ^ A[5] ^ A[10] ^ A[15] ^ A[20];
	C[1] = A[1] ^ A[6] ^ A[11] ^ A[16] ^ A[21];
	C[2] = A[2] ^ A[7] ^ A[12] ^ A[17] ^ A[22];
	C[3] = A[3] ^ A[8] ^ A[13] ^ A[18] ^ A[23];
	C[4] = A[4] ^ A[9] ^ A[14] ^ A[19] ^ A[24];

	/* Compute theta for each column */
	D[0] = C[4] ^ rotl_lane(C[1], 1);
	D[1] = C[0] ^ rotl_lane(C[2], 1);
	D[2] = C[1] ^ rotl_lane(C[3], 1);
	D[3] = C[2] ^ rotl_lane(C[4], 1);
	D[4] = C[3] ^ rotl_lane(C[0], 1);

	/* 1st plane */

	/* Apply theta-rho-pi */
	T[0] = A[0] ^ D[0];
	T[1] = rotl_lane(A[6]  ^ D[1], 44);
	T[2] = rotl_lane(A[12] ^ D[2], 43);
	T[3] = rotl_lane(A[18] ^ D[3], 21);
	T[4] = rotl_lane(A[24] ^ D[4], 14);

	/* Apply chi */
	/* Also apply iota since we are here */
	N[0] = T[0] ^ (~T[1] & T[2]) ^ round_constants[r_idx];
	N[1] = T[1] ^ (~T[2] & T[3]);
	N[2] = T[2] ^ (~T[3] & T[4]);
	N[3] = T[3] ^ (~T[4] & T[0]);
	N[4] = T[4] ^ (~T[0] & T[1]);

	/* 2nd plane */

	T[0] = rotl_lane(A[3]  ^ D[3], 28);
	T[1] = rotl_lane(A[9]  ^ D[4], 20);
	T[2] = rotl_lane(A[10] ^ D[0], 3);
	T[3] = rotl_lane(A[16] ^ D[1], 45);
	T[4] = rotl_lane(A[22] ^ D[2], 61);

	N[5] = T[0] ^ (~T[1] & T[2]);
	N[6] = T[1] ^ (~T[2] & T[3]);
	N[7] = T[2] ^ (~T[3] & T[4]);
	N[8] = T[3] ^ (~T[4] & T[0]);
	N[9] = T[4] ^ (~T[0] & T[1]);

	/* 3rd plane */

	T[0] = rotl_lane(A[1]  ^ D[1], 1);
	T[1] = rotl_lane(A[7]  ^ D[2], 6);
	T[2] = rotl_lane(A[13] ^ D[3], 25);
	T[3] = rotl_lane(A[19] ^ D[4], 8);
	T[4] = rotl_lane(A[20] ^ D[0], 18);

	N[10] = T[0] ^ (~T[1] & T[2]);
	N[11] = T[1] ^ (~T[2] & T[3]);
	N[12] = T[2] ^ (~T[3] & T[4]);
	N[13] = T[3] ^ (~T[4] & T[0]);
	N[14] = T[4] ^ (~T[0] & T[1]);

	/* 4th plane */

	T[0] = rotl_lane(A[4]  ^ D[4], 27);
	T[1] = rotl_lane(A[5]  ^ D[0], 36);
	T[2] = rotl_lane(A[11] ^ D[1], 10);
	T[3] = rotl_lane(A[17] ^ D[2], 15);
	T[4] = rotl_lane(A[23] ^ D[3], 56);

	N[15] = T[0] ^ (~T[1] & T[2]);
	N[16] = T[1] ^ (~T[2] & T[3]);
	N[17] = T[2] ^ (~T[3] & T[4]);
	N[18] = T[3] ^ (~T[4] & T[0]);
	N[19] = T[4] ^ (~T[0] & T[1]);

	/* 5th plane */

	T[0] = rotl_lane(A[2]  ^ D[2], 62);
	T[1] = rotl_lane(A[8]  ^ D[3], 55);
	T[2] = rotl_lane(A[14] ^ D[4], 39);
	T[3] = rotl_lane(A[15] ^ D[0], 41);
	T[4] = rotl_lane(A[21] ^ D[1], 2);

	N[20] = T[0] ^ (~T[1] & T[2]);
	N[21] = T[1] ^ (~T[2] & T[3]);
	N[22] = T[2] ^ (~T[3] & T[4]);
	N[23] = T[3] ^ (~T[4] & T[0]);
	N[24] = T[4] ^ (~T[0] & T[1]);
}

void
keccakf1600_state_permute_intermediateur(k1600_state_t *st)
{
	lane_t N[KECCAK_NUM_LANES] = { 0 };
	int i = 0;

	for (i = 0; i < KECCAK1600_NUM_ROUNDS; i += 2) {
		keccakf1600_round_intermediate_unrolled(st->A, N, i);
		keccakf1600_round_intermediate_unrolled(N, st->A, i + 1);
	}
}
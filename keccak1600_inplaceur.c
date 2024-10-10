// SPDX-License-Identifier: BSD-2-Clause
/*
 * Keccak-f[1600] C Implementation - State permutation
 * Copyright (C) 2019 Nick Kossifidis <mick@ics.forth.gr>
 */

#include "keccak1600.h"

/*
 * This is an unrolled version of the original implementation
 * with some further cleanups to see how well it competes
 * with the compiler's optimizations. It also provides a
 * more compact view of what's going on.
 *
 * When Zbb is enabled, it's the most performant implementation
 * in my tests.
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
keccakf1600_round_inplace_unrolled(lane_t * A, int r_idx)
{
	lane_t C[5] = { 0 };
	lane_t T = 0;

	/* Compute the parity of the columns */
	C[0] = A[0] ^ A[5] ^ A[10] ^ A[15] ^ A[20];
	C[1] = A[1] ^ A[6] ^ A[11] ^ A[16] ^ A[21];
	C[2] = A[2] ^ A[7] ^ A[12] ^ A[17] ^ A[22];
	C[3] = A[3] ^ A[8] ^ A[13] ^ A[18] ^ A[23];
	C[4] = A[4] ^ A[9] ^ A[14] ^ A[19] ^ A[24];

	/* Compute D for each slice (reuse C) */
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
keccakf1600_state_permute_inplaceur(k1600_state_t *st)
{
	int i = 0;
	for (i = 0; i < KECCAK1600_NUM_ROUNDS; i++)
		keccakf1600_round_inplace_unrolled(st->A, i);
}
// SPDX-License-Identifier: BSD-2-Clause
/*
 * Keccak-f[1600] C Implementation - State permutation
 * Copyright (C) 2024 Nick Kossifidis <mick@ics.forth.gr>
 */

#include "keccak1600.h"

/*
 * Same as intermediate_unrolled but using early parity, an optimization
 * comming from section 2.4.1 of "Keccak implementation overview".
 *
 * Instead of computing C[x] for all columns at the begining
 * of each round, do it as soon as columns are available to further improve
 * data locality. The first round is as usual and the last round skips this
 * step. To avoid checking for r_idx on every plane on each round, we have a
 * separate function for the last round.
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
keccakf1600_round_intermediate_unrolled_ep(lane_t *A, lane_t *N, lane_t *C, int r_idx)
{
	lane_t D[5] = { 0 };
	lane_t T[5] = { 0 };

	/* Compute theta for each column
	 * (C already computed)*/
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

	/* Update C */
	C[0] = N[0];
	C[1] = N[1];
	C[2] = N[2];
	C[3] = N[3];
	C[4] = N[4];

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

	C[0] ^= N[5];
	C[1] ^= N[6];
	C[2] ^= N[7];
	C[3] ^= N[8];
	C[4] ^= N[9];

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

	C[0] ^= N[10];
	C[1] ^= N[11];
	C[2] ^= N[12];
	C[3] ^= N[13];
	C[4] ^= N[14];

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

	C[0] ^= N[15];
	C[1] ^= N[16];
	C[2] ^= N[17];
	C[3] ^= N[18];
	C[4] ^= N[19];

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

	C[0] ^= N[20];
	C[1] ^= N[21];
	C[2] ^= N[22];
	C[3] ^= N[23];
	C[4] ^= N[24];
}

static inline void
keccakf1600_round_intermediate_unrolled_ep_last(lane_t *A, lane_t *N, lane_t *C, int r_idx)
{
	lane_t D[5] = { 0 };
	lane_t T[5] = { 0 };

	/* Compute theta for each column
	 * (C already computed)*/
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
keccakf1600_state_permute_intermediateur_ep(k1600_state_t *st)
{
	lane_t N[KECCAK_NUM_LANES] = { 0 };
	lane_t C[5] = { 0 };
	int i = 0;

	C[0] = st->A[0] ^ st->A[5] ^ st->A[10] ^ st->A[15] ^ st->A[20];
	C[1] = st->A[1] ^ st->A[6] ^ st->A[11] ^ st->A[16] ^ st->A[21];
	C[2] = st->A[2] ^ st->A[7] ^ st->A[12] ^ st->A[17] ^ st->A[22];
	C[3] = st->A[3] ^ st->A[8] ^ st->A[13] ^ st->A[18] ^ st->A[23];
	C[4] = st->A[4] ^ st->A[9] ^ st->A[14] ^ st->A[19] ^ st->A[24];

	for (i = 0; i < KECCAK1600_NUM_ROUNDS - 2; i += 2) {
		keccakf1600_round_intermediate_unrolled_ep(st->A, N, C, i);
		keccakf1600_round_intermediate_unrolled_ep(N, st->A, C, i + 1);
	}

	keccakf1600_round_intermediate_unrolled_ep(st->A, N, C, 22);
	keccakf1600_round_intermediate_unrolled_ep_last(N, st->A, C, 23);
}
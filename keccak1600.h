// SPDX-License-Identifier: BSD-2-Clause
/*
 * SHA3 (Keccak-f[1600]) C / RV64 Implementation
 * Copyright (C) 2019 Nick Kossifidis <mick@ics.forth.gr>
 */

/*
 * Keccak's state is a cuboid with 5bits width (x),
 * 5bits height (y) and <lane> bits depth, where a
 * lane is meant to be as large as the CPU's word
 * size. The SHA-3 specification (FIPS 202) however
 * is only defined for keccak-f[1600] which means
 * the state size is fixed to 1600bits and a lane
 * size of 64bits (5*5*64). So it's more fit for
 * 64bit CPUs.
 *
 * For more info check out:
 * https://keccak.team/keccak.html
 */

#include <stdint.h>	/* For typed integers */
#include <stddef.h>	/* For size_t */

#define KECCAK_NUM_COLS		5
#define KECCAK_NUM_ROWS		5
#define KECCAK_NUM_LANES	25	/* KECCAK_NUM_COLS * KECCAK_NUM_ROWS */

typedef uint64_t lane_t;
#define KECCAK1600_LANE_BITS	64
#define	KECCAK1600_LANE_BYTES	8
#define KECCAK1600_NUM_ROUNDS	24	/* 12 + 2^log2(KECCAK1600_LANE_BITS) */
#define KECCAK1600_STATE_SIZE	200	/* 1600 bits (lane_size * num_lanes) */

typedef union {
	/* Rows are allocated from bottom to top
	 * so we get the first row for y = 0 with
	 * x = 0 - 4, then for y = 1 etc. The mapping
	 * from (x,y) to the index of this array is
	 * Index(x, y) = x + 5*y where x and y are
	 * from 0 - 4 (so we % 5 each coordinate
	 * to be sure they stay within range) */
	lane_t A[KECCAK_NUM_LANES];
	uint8_t A_bytes[KECCAK1600_STATE_SIZE];
} k1600_state_t;

/* Left-rotate a lane, keep it like this so that the
 * compiler recognizes it and optimizes it using arch-specific
 * bit manip. instructions. */
static inline lane_t rotl_lane(lane_t val, int times)
{
	return (((val) << (times)) |
		((val) >> (KECCAK1600_LANE_BITS - (times))));
}

/* Used for handling multiple underlying implementations */
typedef void (*keccak1600_spf) (k1600_state_t * st);

void keccakf1600_set_permutation_function(keccak1600_spf func, int lc);
void keccakf1600_oneshot(const void *msg, size_t msg_len, void *md,
			 size_t md_len, uint8_t delim_suffix);

/* Available implementations */
void keccakf1600_state_permute_simple(k1600_state_t *st);
void keccakf1600_state_permute_inplaceur(k1600_state_t *st);
void keccakf1600_state_permute_intermediateur(k1600_state_t *st);
void keccakf1600_state_permute_intermediateur_ep(k1600_state_t *st);
void keccakf1600_state_permute_intermediateur_lc(k1600_state_t *st);

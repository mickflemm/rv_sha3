// SPDX-License-Identifier: BSD-2-Clause
/*
 * Keccak-f[1600] C Implementation - Sponge functions
 * Copyright (C) 2019 Nick Kossifidis <mick@ics.forth.gr>
 */

#include <stddef.h>		/* For size_t */
#include <string.h>		/* For memcpy() */
#include "keccak1600.h"

/******************\
* SPONGE FUNCTIONS *
\******************/

static inline void
keccakf1600_absorb_lanes(k1600_state_t * st, const void *msg,
			 int num_lanes, int *msg_off)
{
	int i = 0;

	for (i = 0; i < num_lanes; i++, (*msg_off) += KECCAK1600_LANE_BYTES)
		st->A[i] ^= (*((lane_t *) (msg + (*msg_off))));

	keccakf1600_state_permute(st);
}

static void
keccakf1600_absorb(k1600_state_t * st, const void *msg, size_t msg_len,
		   size_t md_len, uint8_t delim_suffix)
{
	int capacity_bytes = 2 * md_len;
	int rate_bytes = KECCAK1600_STATE_SIZE - capacity_bytes;
	int num_blocks = msg_len / rate_bytes;
	int lanes_per_block = rate_bytes / KECCAK1600_LANE_BYTES;
	int msg_off = 0;
	int block_off = 0;
	int i = 0;

	/* Absorb message */

	/* Blocks are multiples of lane size
	 * so absorb a lane at a time (instead of
	 * a byte at a time) to speed things up */
	for (i = 0; i < num_blocks; i++)
		keccakf1600_absorb_lanes(st, msg, lanes_per_block, &msg_off);

	/* Handle any remaining bytes */
	for (i = msg_off; i < msg_len; i++) {
		st->A_bytes[block_off++] ^= ((const uint8_t *)msg)[i];
		if (block_off == rate_bytes) {
			keccakf1600_state_permute(st);
			block_off = 0;
		}
	}

	/* Absorb padding */
	/* For delim_suffix check out
	 * https://keccak.team/keccak_specs_summary.html */
	st->A_bytes[block_off] ^= delim_suffix;

	/* The delimiter is at the end of the block, we need
	 * another block for the second bit of padding, absorb
	 * this one and work on the next */
	if ((delim_suffix & 0x80) && (block_off == (rate_bytes - 1)))
		keccakf1600_state_permute(st);

	st->A_bytes[rate_bytes - 1] ^= 0x80;
	keccakf1600_state_permute(st);
}

static void
keccakf1600_squeeze(k1600_state_t * st, void *md, size_t md_len)
{
	int capacity_bytes = 2 * md_len;
	int rate_bytes = KECCAK1600_STATE_SIZE - capacity_bytes;
	int block_len = 0;
	char *md_off = md;
	int i = md_len;
	int c = 0;

	while (i > 0) {
		block_len = (i < rate_bytes) ? i : rate_bytes;
		memcpy(md_off, st->A_bytes, block_len);

		md_off += block_len;
		i -= block_len;

		/* Squeeze another block out of the state */
		if (i > 0)
			keccakf1600_state_permute(st);
	}
}


/*************\
* ENTRY POINT *
\*************/

void
keccakf1600_oneshot(const void *msg, size_t msg_len, void *md,
		    size_t md_len, uint8_t delim_suffix)
{
	k1600_state_t st = { 0 };
	keccakf1600_absorb(&st, msg, msg_len, md_len, delim_suffix);
	keccakf1600_squeeze(&st, md, md_len);
}

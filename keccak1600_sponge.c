// SPDX-License-Identifier: BSD-2-Clause
/*
 * Keccak-f[1600] C Implementation - Sponge functions
 * Copyright (C) 2019 Nick Kossifidis <mick@ics.forth.gr>
 */

#include "keccak1600.h"
#include <string.h>		/* For memcpy() */
#include <stdbool.h>		/* For bool */

/*******************************\
* WRAPPER FOR STATE PERMUTATION *
\*******************************/
/* Stub function so that we don't check for NULL every time */
static void stub_state_permute(k1600_state_t * st)  { return; }
static keccak1600_spf keccakf1600_state_permute = &stub_state_permute;
static bool use_lc = false;

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
keccakf1600_squeeze(k1600_state_t *st, void *md, size_t md_len)
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

		if (use_lc) {
			/* Apply P mask to output */
			int lanes_out = block_len / 8;
			lane_t *md_lanes = (lane_t *) md_off;
			if (lanes_out > 1)
				md_lanes[1] = ~md_lanes[1];
			if (lanes_out > 2)
				md_lanes[2] = ~md_lanes[2];
			if (lanes_out > 8)
				md_lanes[8] = ~md_lanes[8];
			if (lanes_out > 12)
				md_lanes[12] = ~md_lanes[12];
			if (lanes_out > 17)
				md_lanes[17] = ~md_lanes[17];
			if (lanes_out > 20)
				md_lanes[20] = ~md_lanes[20];
		}
		
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
keccakf1600_set_permutation_function(keccak1600_spf func, int lc)
{
	keccakf1600_state_permute = func;
	if (lc)
		use_lc = true;
	else
		use_lc = false;
}

void
keccakf1600_oneshot(const void *msg, size_t msg_len, void *md,
		    size_t md_len, uint8_t delim_suffix)
{
	k1600_state_t st = { 0 };

	/* When doing lane complementing, operate on a
	 * partialy inverted state. */
	if (use_lc) {
		st.A[1] = ~0ULL;
		st.A[2] = ~0ULL;
		st.A[8] = ~0ULL;
		st.A[12] = ~0ULL;
		st.A[17] = ~0ULL;
		st.A[20] = ~0ULL;
	}
	keccakf1600_absorb(&st, msg, msg_len, md_len, delim_suffix);
	keccakf1600_squeeze(&st, md, md_len);
}

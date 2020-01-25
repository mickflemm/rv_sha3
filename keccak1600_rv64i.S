// SPDX-License-Identifier: BSD-2-Clause
/*
 * Keccak-f[1600] RV64I Implementation - State permutation
 * Copyright (C) 2019 Nick Kossifidis <mick@ics.forth.gr>
 */

#define KECCAK_NUM_ROWS		5
#define KECCAK_NUM_COLS		5
#define KECCAK_NUM_LANES	25

#define KECCAK1600_LANE_BITS	64
#define	KECCAK1600_LANE_BYTES	8
#define KECCAK1600_NUM_ROUNDS	24	/* 12 + 2^log2(KECCAK1600_LANE_BITS) */
#define KECCAK1600_STATE_SIZE	200	/* 1600 bits (lane_size * num_lanes) */

.data

.align 16,0
round_constants:
    .dword 0x0000000000000001
    .dword 0x0000000000008082
    .dword 0x800000000000808a
    .dword 0x8000000080008000
    .dword 0x000000000000808b
    .dword 0x0000000080000001
    .dword 0x8000000080008081
    .dword 0x8000000000008009
    .dword 0x000000000000008a
    .dword 0x0000000000000088
    .dword 0x0000000080008009
    .dword 0x000000008000000a
    .dword 0x000000008000808b
    .dword 0x800000000000008b
    .dword 0x8000000000008089
    .dword 0x8000000000008003
    .dword 0x8000000000008002
    .dword 0x8000000000000080
    .dword 0x000000000000800a
    .dword 0x800000008000000a
    .dword 0x8000000080008081
    .dword 0x8000000000008080
    .dword 0x0000000080000001
    .dword 0x8000000080008008

.text

/*
 * Reserved registers:
 *
 * s11 -> Return address
 * s10 -> Round counter
 * s9  -> KECCAK1600_NUM_ROUNDS
 * s8  -> KECCAK1600_LANE_BITS
 * s7  -> Pointer to round_constants
 */


/*********\
* HELPERS *
\*********/

/*
 * Left-rotate a lane
 * s8  -> KECCAK1600_LANE_BITS
 */
.macro ROTL_LANE _out_lane _in_lane _times
	li	t1, \_times
	sll	t0, \_in_lane, t1
	sub	t1, s8, t1
	srl	t1, \_in_lane, t1
	or	t0, t0, t1
	mv	\_out_lane, t0
.endm

/*
 * Load 5 lanes from A to a3 - a7
 * a0 -> Pointer to A
 */
.macro LOAD_5LANES _a _b _c _d _e
	ld	a3, (8 * \_a)(a0)
	ld	a4, (8 * \_b)(a0)
	ld	a5, (8 * \_c)(a0)
	ld	a6, (8 * \_d)(a0)
	ld	a7, (8 * \_e)(a0)
.endm

/*
 * Store 5 lanes from a3 - a7 to N
 * a1 -> Pointer to N
 */
.macro STORE_5LANES _a _b _c _d _e
	sd	a3, (8 * \_a)(a1)
	sd	a4, (8 * \_b)(a1)
	sd	a5, (8 * \_c)(a1)
	sd	a6, (8 * \_d)(a1)
	sd	a7, (8 * \_e)(a1)
.endm


/**********************\
* KECCAK STEP MAPPINGS *
\**********************/

/*
 * A combined theta_rho step
 * out_lane = rotl64(in_lane ^ D_val, rho_off)
 * in_lane comes from the loaded 
 */
.macro THETA_RHO_STEP _out_lane, _in_lane, _d_val, _rho_off
	xor		\_out_lane, \_in_lane, \_d_val
	ROTL_LANE 	\_out_lane, \_out_lane, \_rho_off
.endm

/*
 * A Chi step on the plane
 *
 * Reuse a3 - a7 so that we can
 * do iota on the 1st lane of
 * the 1st plane and also write the
 * plane in one take
 *
 * a3 - a6 -> N[0]
 * t2 - t6 -> T[0]
 */
.macro CHI
	not	a3, t3
	and	a3, a3, t4
	xor	a3, a3, t2

	not	a4, t4
	and	a4, a4, t5
	xor	a4, a4, t3

	not	a5, t5
	and	a5, a5, t6
	xor	a5, a5, t4

	not	a6, t6
	and	a6, a6, t2
	xor	a6, a6, t5

	not	a7, t2
	and	a7, a7, t3
	xor	a7, a7, t6
.endm

/*
 * Apply iota by xoring the
 * first lane of the first
 * plane, stored in a3 after
 * CHI, with the round constant
 * for this round
 *
 * s10 -> Round counter
 * s7 -> Pointer to round_constants
 */
.macro IOTA
	slli	t0, s10, 3
	add	t0, t0, s7
	ld	t1, 0(t0)
	xor	a3, a3, t1
.endm

/*
 * This is the equivalent of keccakf1600_round_unrolled_cf
 * since we don't have enough registers to play with.
 *
 * a0 -> Pointer to state A
 * a1 -> Pointer to state N
 * a2 -> Round index
 * s1 - s5 -> C[0] / D[0]
 * t2 - t6 -> T[0]
 */
.align 8,0
.func round_unrolled_cf
.global round_unrolled_cf
round_unrolled_cf:

	/* Compute parity of columns and place the
	 * results on s1 - s5 */

	ld 	t0, 0(a0)	/* C[0] = A[0] */
	ld 	t1, 40(a0)
	ld	t2, 80(a0)
	ld	t3, 120(a0)
	ld	t4, 160(a0)
	xor 	t0, t0, t1	/* ^ A[5] */
	xor	t0, t0, t2	/* ^ A[10] */
	xor	t0, t0, t3	/* ^ A[15] */
	xor	t0, t0, t4	/* ^ A[20] */
	mv	s1, t0

	ld 	t0, 8(a0)	/* C[1] = A[1] */
	ld 	t1, 48(a0)
	ld	t2, 88(a0)
	ld	t3, 128(a0)
	ld	t4, 168(a0)
	xor 	t0, t0, t1	/* ^ A[6] */
	xor	t0, t0, t2	/* ^ A[11] */
	xor	t0, t0, t3	/* ^ A[16] */
	xor	t0, t0, t4	/* ^ A[21] */
	mv	s2, t0

	ld 	t0, 16(a0)	/* C[2] = A[2] */
	ld 	t1, 56(a0)
	ld	t2, 96(a0)
	ld	t3, 136(a0)
	ld	t4, 176(a0)
	xor 	t0, t0, t1	/* ^ A[7] */
	xor	t0, t0, t2	/* ^ A[12] */
	xor	t0, t0, t3	/* ^ A[17] */
	xor	t0, t0, t4	/* ^ A[22] */
	mv	s3, t0

	ld 	t0, 24(a0)	/* C[3] = A[3] */
	ld 	t1, 64(a0)
	ld	t2, 104(a0)
	ld	t3, 144(a0)
	ld	t4, 184(a0)
	xor 	t0, t0, t1	/* ^ A[8] */
	xor	t0, t0, t2	/* ^ A[13] */
	xor	t0, t0, t3	/* ^ A[18] */
	xor	t0, t0, t4	/* ^ A[23] */
	mv	s4, t0

	ld 	t0, 32(a0)	/* C[4] = A[4] */
	ld 	t1, 72(a0)
	ld	t2, 112(a0)
	ld	t3, 152(a0)
	ld	t4, 192(a0)
	xor 	t0, t0, t1	/* ^ A[9] */
	xor	t0, t0, t2	/* ^ A[14] */
	xor	t0, t0, t3	/* ^ A[19] */
	xor	t0, t0, t4	/* ^ A[24] */
	mv	s5, t0

	/* Calculate theta and overwrite s1 - s5 with
	 * the D values */

	/* T[0] = C[4] */
	mv	t2, s5

	/* C[4] ^= rotl_lane(C[1], 1);	-> D[0] */
	ROTL_LANE	t3, s2, 1
	xor		s5, s5, t3

	/* C[1] ^= rotl_lane(C[3], 1);	-> D[2] */
	ROTL_LANE	t3, s4, 1
	xor		s2, s2, t3

	/* C[3] ^= rotl_lane(C[0], 1);	-> D[4] */
	ROTL_LANE	t3, s1, 1
	xor		s4, s4, t3

	/* C[0] ^= rotl_lane(C[2], 1);	-> D[1] */
	ROTL_LANE	t3, s3, 1
	xor		s1, s1, t3

	/* C[2] ^= rotl_lane(T[0], 1);	-> D[3] */
	ROTL_LANE	t3, t2, 1
	xor		s3, s3, t3


	/* Plane by plane processing,
	 * T is stored on t2 - t6 */

	LOAD_5LANES	0, 6, 12, 18, 24

	xor		t2, a3, s5
	THETA_RHO_STEP	t3, a4, s1, 44
	THETA_RHO_STEP	t4, a5, s2, 43
	THETA_RHO_STEP	t5, a6, s3, 21
	THETA_RHO_STEP	t6, a7, s4, 14

	/* N will be at a3 - a7 */
	CHI

	/* This is the first plane, apply iota
	 * on N[0] (a3) since we are here */
	IOTA

	STORE_5LANES	0, 1, 2, 3, 4

	/* 2nd plane */
	LOAD_5LANES	3, 9, 10, 16, 22
	THETA_RHO_STEP	t2, a3, s3, 28
	THETA_RHO_STEP	t3, a4, s4, 20
	THETA_RHO_STEP	t4, a5, s5, 3
	THETA_RHO_STEP	t5, a6, s1, 45
	THETA_RHO_STEP	t6, a7, s2, 61
	CHI
	STORE_5LANES	5, 6, 7, 8, 9

	/* 3rd plane */
	LOAD_5LANES	1, 7, 13, 19, 20
	THETA_RHO_STEP	t2, a3, s1, 1
	THETA_RHO_STEP	t3, a4, s2, 6
	THETA_RHO_STEP	t4, a5, s3, 25
	THETA_RHO_STEP	t5, a6, s4, 8
	THETA_RHO_STEP	t6, a7, s5, 18
	CHI
	STORE_5LANES	10, 11, 12, 13, 14

	/* 4rth plane */
	LOAD_5LANES	4, 5, 11, 17, 23
	THETA_RHO_STEP	t2, a3, s4, 27
	THETA_RHO_STEP	t3, a4, s5, 36
	THETA_RHO_STEP	t4, a5, s1, 10
	THETA_RHO_STEP	t5, a6, s2, 15
	THETA_RHO_STEP	t6, a7, s3, 56
	CHI
	STORE_5LANES	15, 16, 17, 18, 19

	/* 5th plane */
	LOAD_5LANES	2, 8, 14, 15, 21
	THETA_RHO_STEP	t2, a3, s2, 62
	THETA_RHO_STEP	t3, a4, s3, 55
	THETA_RHO_STEP	t4, a5, s4, 39
	THETA_RHO_STEP	t5, a6, s5, 41
	THETA_RHO_STEP	t6, a7, s1, 2
	CHI
	STORE_5LANES	20, 21, 22, 23, 24

	ret
.endfunc


/*************\
* ENTRY POINT *
\*************/

/*
 * a0 -> Pointer to A
 */
.align 8,0
.func keccakf1600_state_permute
.global keccakf1600_state_permute
keccakf1600_state_permute:
	/* Save s registers on the stack and make
	 * room for the intermediate state T */
	addi	sp, sp, -304
	sd	s11, 88(sp)
	sd	s10, 80(sp)
	sd	s9, 72(sp)
	sd	s8, 64(sp)
	sd	s7, 56(sp)
	sd	s6, 48(sp)
	sd	s5, 40(sp)
	sd	s4, 32(sp)
	sd	s3, 24(sp)
	sd	s2, 16(sp)
	sd	s1, 8(sp)
	sd	s0, 0(sp)

	/* Update frame pointer */
	addi	s0, sp, 304

	/* Save return address since
	 * we'll be making calls from
	 * here as well */
	mv 	s11, ra

	/* Initialize state */
	mv	s10, zero
	li	s9, KECCAK1600_NUM_ROUNDS
	li	s8, KECCAK1600_LANE_BITS
	la	s7, round_constants
	addi	a1, sp, 96

1:
	call	round_unrolled_cf
	mv	t0, a0
	mv	a0, a1
	mv	a1, t0
	addi	s10, s10, 1
	blt	s10, s9, 1b

	/* Restore return address and
	 * s registers from stack */
	mv 	ra, s11
	ld	s0, 0(sp)
	ld	s1, 8(sp)
	ld	s2, 16(sp)
	ld	s3, 24(sp)
	ld	s4, 32(sp)
	ld	s5, 40(sp)
	ld	s6, 48(sp)
	ld	s7, 56(sp)
	ld	s8, 64(sp)
	ld	s9, 72(sp)
	ld	s10, 80(sp)
	ld	s11, 88(sp)
	addi	sp, sp, 304
	ret
.endfunc

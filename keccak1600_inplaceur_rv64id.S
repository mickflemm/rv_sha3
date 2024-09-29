// SPDX-License-Identifier: BSD-2-Clause
/*
 * Keccak-f[1600] RV64ID Implementation - State permutation
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

.align 3
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
 * a0  -> Pointer to state
 * s0  -> Frame pointer
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

/* Note that we only need to save / restore
 * FP regs to / from the stack, in case the
 * rest of our code uses them (for example
 * for the sha3_test), if not we can skip this
 * part.
 *
 * Also note that we don't need to save/restore
 * all of them according to the ABI, so this approach
 * uses less memory than intermediateur, even when we
 * do save/restore FP regs.
 */
.macro SAVE_FP_REGS _sp_off
	fsd	fa7, (\_sp_off + 128)(sp)
	fsd	fa6, (\_sp_off + 120)(sp)
	fsd	fa5, (\_sp_off + 112)(sp)
	fsd	fa4, (\_sp_off + 104)(sp)
	fsd	fa3, (\_sp_off + 96)(sp)
	fsd	fa2, (\_sp_off + 88)(sp)
	fsd	fa1, (\_sp_off + 80)(sp)
	fsd	fa0, (\_sp_off + 72)(sp)
	fsd	fs8, (\_sp_off + 64)(sp)
	fsd	fs7, (\_sp_off + 56)(sp)
	fsd	fs6, (\_sp_off + 48)(sp)
	fsd	fs5, (\_sp_off + 40)(sp)
	fsd	fs4, (\_sp_off + 32)(sp)
	fsd	fs3, (\_sp_off + 24)(sp)
	fsd	fs2, (\_sp_off + 16)(sp)
	fsd	fs1, (\_sp_off + 8)(sp)
	fsd	fs0, (\_sp_off)(sp)
.endm

.macro RESTORE_FP_REGS _sp_off
	fld	fs0, (\_sp_off)(sp)
	fld	fs1, (\_sp_off + 8)(sp)
	fld	fs2, (\_sp_off + 16)(sp)
	fld	fs3, (\_sp_off + 24)(sp)
	fld	fs4, (\_sp_off + 32)(sp)
	fld	fs5, (\_sp_off + 40)(sp)
	fld	fs6, (\_sp_off + 48)(sp)
	fld	fs7, (\_sp_off + 56)(sp)
	fld	fs8, (\_sp_off + 64)(sp)
	fld	fa0, (\_sp_off + 72)(sp)
	fld	fa1, (\_sp_off + 80)(sp)
	fld	fa2, (\_sp_off + 88)(sp)
	fld	fa3, (\_sp_off + 96)(sp)
	fld	fa4, (\_sp_off + 104)(sp)
	fld	fa5, (\_sp_off + 112)(sp)
	fld	fa6, (\_sp_off + 120)(sp)
	fld	fa7, (\_sp_off + 128)(sp)
.endm

.macro GET_LANE _out, _num
	fmv.x.d	\_out, f\_num
.endm

.macro SET_LANE _num, _in
	fmv.d.x	f\_num, \_in
.endm

.macro LOAD_LANE _num
	fld	f\_num, (8 * \_num)(a0)
.endm

.macro STORE_LANE _num
	fsd	f\_num, (8 * \_num)(a0)
.endm

.macro LOAD_STATE
	LOAD_LANE	0
	LOAD_LANE	1
	LOAD_LANE	2
	LOAD_LANE	3
	LOAD_LANE	4
	LOAD_LANE	5
	LOAD_LANE	6
	LOAD_LANE	7
	LOAD_LANE	8
	LOAD_LANE	9
	LOAD_LANE	10
	LOAD_LANE	11
	LOAD_LANE	12
	LOAD_LANE	13
	LOAD_LANE	14
	LOAD_LANE	15
	LOAD_LANE	16
	LOAD_LANE	17
	LOAD_LANE	18
	LOAD_LANE	19
	LOAD_LANE	20
	LOAD_LANE	21
	LOAD_LANE	22
	LOAD_LANE	23
	LOAD_LANE	24
.endm

.macro STORE_STATE
	STORE_LANE	0
	STORE_LANE	1
	STORE_LANE	2
	STORE_LANE	3
	STORE_LANE	4
	STORE_LANE	5
	STORE_LANE	6
	STORE_LANE	7
	STORE_LANE	8
	STORE_LANE	9
	STORE_LANE	10
	STORE_LANE	11
	STORE_LANE	12
	STORE_LANE	13
	STORE_LANE	14
	STORE_LANE	15
	STORE_LANE	16
	STORE_LANE	17
	STORE_LANE	18
	STORE_LANE	19
	STORE_LANE	20
	STORE_LANE	21
	STORE_LANE	22
	STORE_LANE	23
	STORE_LANE	24
.endm

/*
 * Calculate parity of column _col and store value to _out
 */
.macro COLUMN_PARITY _out, _a, _b, _c, _d, _e
	GET_LANE t0, \_a
	GET_LANE t1, \_b
	GET_LANE t2, \_c
	GET_LANE t3, \_d
	GET_LANE t4, \_e

	xor	t1, t1, t0
	xor	t1, t1, t2
	xor	t1, t1, t3
	xor	t1, t1, t4
	mv	\_out, t1
.endm

/**********************\
* KECCAK STEP MAPPINGS *
\**********************/

/*
 * A combined theta_rho step
 * out_lane = rotl64(in_lane ^ D_val, rho_off)
 * both in_lane and out_lane are fp regs
 */
.macro THETA_RHO_STEP _out_lane, _in_lane, _d_val, _rho_off
	GET_LANE 	t3, \_in_lane
	xor		t3, t3, \_d_val
	ROTL_LANE 	t3, t3, \_rho_off
	SET_LANE	\_out_lane, t3
.endm

/*
 * A Chi step on the plane
 */
.macro CHI _a _b _c _d _e
	GET_LANE t0, \_a
	GET_LANE t1, \_b
	GET_LANE t2, \_c
	GET_LANE t3, \_d
	GET_LANE t4, \_e

	not	t5, t1
	and	t5, t5, t2
	xor	a1, t5, t0

	not	t5, t2
	and	t5, t5, t3
	xor	a2, t5, t1

	not	t5, t3
	and	t5, t5, t4
	xor	a3, t5, t2

	not	t5, t4
	and	t5, t5, t0
	xor	a4, t5, t3

	not	t5, t0
	and	t5, t5, t1
	xor	a5, t5, t4

	SET_LANE \_a, a1
	SET_LANE \_b, a2
	SET_LANE \_c, a3
	SET_LANE \_d, a4
	SET_LANE \_e, a5
.endm

/*
 * Apply iota by xoring the
 * first lane of the first
 * plane, with the round
 * constant for this round
 *
 * s10 -> Round counter
 * s7 -> Pointer to round_constants
 */
.macro IOTA
	slli	t0, s10, 3
	add	t0, t0, s7
	ld	t1, 0(t0)
	GET_LANE t0, 0
	xor	t0, t0, t1
	SET_LANE 0, t0
.endm

/*************\
* ENTRY POINT *
\*************/

/*
 * This is the equivalent of keccakf1600_inplaceur
 * using fp registers to hold the state, reducing
 * memory accesses. This way we only read/write state
 * from/to memory once, and only use registers during
 * state permutation.
 *
 * a0 -> Pointer to A
 * f0 - f24 -> State
 * t2,f28 -> temp
 * s1 - s5 -> C[] / D[]
 */
.align 3
.func keccakf1600_state_permute_inplaceur_rv64id
.global keccakf1600_state_permute_inplaceur_rv64id
keccakf1600_state_permute_inplaceur_rv64id:
	/* Save s registers on the stack and make
	 * room for the intermediate state T */
	addi	sp, sp, -232
	SAVE_FP_REGS 96
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
	addi	s0, sp, 232

	/* Save return address since
	 * we'll be making calls from
	 * here as well */
	mv 	s11, ra

	/* Initialize state */
	mv	s10, zero
	li	s9, KECCAK1600_NUM_ROUNDS
	li	s8, KECCAK1600_LANE_BITS
	la	s7, round_constants

	LOAD_STATE
1:
	/* Compute parity of columns and place the
	 * results on s1 - s5 */
	/* C[i] = A[i] ^ A[i + 5] ^ A[i + 10] ^ A[i + 15] ^ A[i + 20] */
	COLUMN_PARITY s1, 0, 5, 10, 15, 20
	COLUMN_PARITY s2, 1, 6, 11, 16, 21
	COLUMN_PARITY s3, 2, 7, 12, 17, 22
	COLUMN_PARITY s4, 3, 8, 13, 18, 23
	COLUMN_PARITY s5, 4, 9, 14, 19, 24

	/* Calculate theta and overwrite s1 - s5 with
	 * the D values */

	/* T = C[4] */
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

	/* C[2] ^= rotl_lane(C[4], 1);	-> D[3] */
	ROTL_LANE	t3, t2, 1
	xor		s3, s3, t3

	/* Apply theta-rho-pli in-place */

	/* A[0] ^= D[0] */
	GET_LANE t0, 0
	xor	t0, t0, s5
	SET_LANE 0, t0

	/* T = A[1] */
	GET_LANE t2, 1
	THETA_RHO_STEP 1, 6, s1, 44
	THETA_RHO_STEP 6, 9, s4, 20
	THETA_RHO_STEP 9, 22, s2, 61
	THETA_RHO_STEP 22, 14, s4, 39
	THETA_RHO_STEP 14, 20, s5, 18
	THETA_RHO_STEP 20, 2, s2, 62
	THETA_RHO_STEP 2, 12, s2, 43
	THETA_RHO_STEP 12, 13, s3, 25
	THETA_RHO_STEP 13, 19, s4, 8
	THETA_RHO_STEP 19, 23, s3, 56
	THETA_RHO_STEP 23, 15, s5, 41
	THETA_RHO_STEP 15, 4, s4, 27
	THETA_RHO_STEP 4, 24, s4, 14
	THETA_RHO_STEP 24, 21, s1, 2
	THETA_RHO_STEP 21, 8, s3, 55
	THETA_RHO_STEP 8, 16, s1, 45
	THETA_RHO_STEP 16, 5, s5, 36
	THETA_RHO_STEP 5, 3, s3, 28
	THETA_RHO_STEP 3, 18, s3, 21
	THETA_RHO_STEP 18, 17, s2, 15
	THETA_RHO_STEP 17, 11, s1, 10
	THETA_RHO_STEP 11, 7, s2, 6
	THETA_RHO_STEP 7, 10, s5, 3
	/* Use f28 for passing T to THETA_RHO_STEP */
	SET_LANE 28, t2
	THETA_RHO_STEP 10, 28, s1, 1

	CHI 0 1 2 3 4
	CHI 5 6 7 8 9
	CHI 10 11 12 13 14
	CHI 15 16 17 18 19
	CHI 20 21 22 23 24

	IOTA

	addi	s10, s10, 1
	blt	s10, s9, 1b
	STORE_STATE

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
	RESTORE_FP_REGS 96
	addi	sp, sp, 232
	ret
.endfunc

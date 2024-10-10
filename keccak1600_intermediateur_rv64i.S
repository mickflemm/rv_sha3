// SPDX-License-Identifier: BSD-2-Clause
/*
 * Keccak-f[1600] RV64I Implementation - State permutation
 * Copyright (C) 2024 Nick Kossifidis <mick@ics.forth.gr>
 */


/*
 * This is the equivalent of keccakf1600_intermediateur
 * since we don't have enough registers to play with.
 * In my tests this is slightly smaller than the compiled
 * C version, uses less memory, and is also a bit faster,
 * but overall the differnce is minimal, I have it here
 * mostly for reference.
 *
 * Allocated registers:
 * a0 -> Pointer to state A
 * sp -> Pointer to state N (since N is allocated on the stack)
 * t1 - t5 -> T[]
 * a1 - a5 -> a[]/A[]/N[]
 * s1 - s5 -> C[]/D[]
 * a6 -> Round counter
 * a7 -> Pointer to round_constants
 */


#define KECCAK1600_NUM_ROUNDS	24

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

/*********\
* HELPERS *
\*********/

.macro _ROTL _out, _in, _times
	#if defined(__riscv_zbb)
	rori	\_out, \_in, (64 - \_times)
	#else
	slli	t0, \_in, \_times
	srli	t6, \_in, (64 - \_times)
	or	\_out, t0, t6
	#endif
.endm

.macro _ANDN _out, _a1, _a2
	#if defined(__riscv_zbb)
	andn	\_out, \_a2, \_a1
	#else
	not	\_out, \_a1
	and	\_out, \_out, \_a2
	#endif
.endm


/*
 * Load 5 lanes from A to a1 - a5
 * a0 -> Pointer to A
 */
.macro LOAD_5LANES _a _b _c _d _e
	ld	a1, (8 * \_a)(a0)
	ld	a2, (8 * \_b)(a0)
	ld	a3, (8 * \_c)(a0)
	ld	a4, (8 * \_d)(a0)
	ld	a5, (8 * \_e)(a0)
.endm

/*
 * Store 5 lanes from a1 - a5 to N
 * sp -> Pointer to N
 */
.macro STORE_5LANES _a _b _c _d _e
	sd	a1, (8 * \_a)(sp)
	sd	a2, (8 * \_b)(sp)
	sd	a3, (8 * \_c)(sp)
	sd	a4, (8 * \_d)(sp)
	sd	a5, (8 * \_e)(sp)
.endm

/*
 * Calculate parity of column _col and store value to _out
 * C[i] = A[i] ^ A[i + 5] ^ A[i + 10] ^ A[i + 15] ^ A[i + 20]
 */
.macro COLUMN_PARITY _out _col
	LOAD_5LANES \_col, (\_col + 5), (\_col + 10), (\_col + 15), (\_col + 20)
	xor	a1, a1, a2
	xor	a2, a3, a4
	xor	a3, a5, a1
	xor	\_out, a3, a2
.endm


/**********************\
* KECCAK STEP MAPPINGS *
\**********************/

/*
 * A combined theta_rho step
 * out_lane = rotl64(in_lane ^ D_val, rho_off)
 * note that in_lane is already pi-mapped
 */
.macro THETA_RHO_STEP _out_lane, _in_lane, _d_val, _rho_off
	xor	\_out_lane, \_in_lane, \_d_val
	_ROTL	\_out_lane, \_out_lane, \_rho_off
.endm

/*
 * A Chi step on the plane
 * A[i + y_offset] ^= (~a[i + 1] & a[i + 2])
 *
 * Here a is on T[] and A[i : i + 4] on a1 - a5.
 *
 * Reuse a1 - a5 so that we store them afterwards
 * to N. Also group _ANDNs and xors together to
 * prevent stalls when possible (e.g. when zbb
 * is available).
 */
.macro CHI
	_ANDN	a1, t2, t3
	_ANDN	a2, t3, t4
	_ANDN	a3, t4, t5
	_ANDN	a4, t5, t1
	_ANDN	a5, t1, t2

	xor	a1, a1, t1
	xor	a2, a2, t2
	xor	a3, a3, t3
	xor	a4, a4, t4
	xor	a5, a5, t5
.endm

/*
 * Apply iota on A[0] (a1) by xoring it
 * with the round constant for this round
 *
 * a6 -> Round counter
 * a7 -> Pointer to round_constants
 */
.macro IOTA
	slli	t0, a6, 3
	add	t0, t0, a7
	ld	t0, 0(t0)
	xor	a1, a1, t0
.endm


/*************\
* ENTRY POINT *
\*************/

.align 3
.func keccakf1600_state_permute_intermediateur_rv64i
.global keccakf1600_state_permute_intermediateur_rv64i
keccakf1600_state_permute_intermediateur_rv64i:
	/* Save s registers on the stack and make
	 * room for the intermediate state N. No
	 * need to save ra since we won't be calling
	 * any functions from here, and also ignore
	 * the frame pointer. */
	addi	sp, sp, -240
	sd	s1, 232(sp)
	sd	s2, 224(sp)
	sd	s3, 216(sp)
	sd	s4, 208(sp)
	sd	s5, 200(sp)

	/* Initialize state */
	mv	a6, zero
	la	a7, round_constants

1:
	/* Compute parity of columns and place the
	 * results on s1 - s5 */
	COLUMN_PARITY s1 0
	COLUMN_PARITY s2 1
	COLUMN_PARITY s3 2
	COLUMN_PARITY s4 3
	COLUMN_PARITY s5 4

	/* Calculate theta and overwrite s1 - s5 with D[]
	 * C[i -1] ^= rotl_lane(C[i + 1], 1) -> D[i] */

	_ROTL	t1, s1, 1
	_ROTL	t2, s2, 1
	_ROTL	t3, s3, 1
	_ROTL	t4, s4, 1
	_ROTL	t5, s5, 1

	xor	s5, s5, t2	// D[0]
	xor	s1, s1, t3	// D[1]
	xor	s2, s2, t4	// D[2]
	xor	s3, s3, t5	// D[3]
	xor	s4, s4, t1	// D[4]

	/* 1st plane */
	LOAD_5LANES	0, 6, 12, 18, 24
	xor		t1, a1, s5
	THETA_RHO_STEP	t2, a2, s1, 44
	THETA_RHO_STEP	t3, a3, s2, 43
	THETA_RHO_STEP	t4, a4, s3, 21
	THETA_RHO_STEP	t5, a5, s4, 14
	CHI
	/* This is the first plane, apply iota
	 * on N[0] (a3) since we are here */
	IOTA
	STORE_5LANES	0, 1, 2, 3, 4

	/* 2nd plane */
	LOAD_5LANES	3, 9, 10, 16, 22
	THETA_RHO_STEP	t1, a1, s3, 28
	THETA_RHO_STEP	t2, a2, s4, 20
	THETA_RHO_STEP	t3, a3, s5, 3
	THETA_RHO_STEP	t4, a4, s1, 45
	THETA_RHO_STEP	t5, a5, s2, 61
	CHI
	STORE_5LANES	5, 6, 7, 8, 9

	/* 3rd plane */
	LOAD_5LANES	1, 7, 13, 19, 20
	THETA_RHO_STEP	t1, a1, s1, 1
	THETA_RHO_STEP	t2, a2, s2, 6
	THETA_RHO_STEP	t3, a3, s3, 25
	THETA_RHO_STEP	t4, a4, s4, 8
	THETA_RHO_STEP	t5, a5, s5, 18
	CHI
	STORE_5LANES	10, 11, 12, 13, 14

	/* 4rth plane */
	LOAD_5LANES	4, 5, 11, 17, 23
	THETA_RHO_STEP	t1, a1, s4, 27
	THETA_RHO_STEP	t2, a2, s5, 36
	THETA_RHO_STEP	t3, a3, s1, 10
	THETA_RHO_STEP	t4, a4, s2, 15
	THETA_RHO_STEP	t5, a5, s3, 56
	CHI
	STORE_5LANES	15, 16, 17, 18, 19

	/* 5th plane */
	LOAD_5LANES	2, 8, 14, 15, 21
	THETA_RHO_STEP	t1, a1, s2, 62
	THETA_RHO_STEP	t2, a2, s3, 55
	THETA_RHO_STEP	t3, a3, s4, 39
	THETA_RHO_STEP	t4, a4, s5, 41
	THETA_RHO_STEP	t5, a5, s1, 2
	CHI
	STORE_5LANES	20, 21, 22, 23, 24

	/* Swap a0 <-> sp, since we have an
	 * even number of rounds, we'll end up
	 * with the correct value on sp. */
	mv	t0, a0
	mv	a0, sp
	mv	sp, t0
	addi	a6, a6, 1
	li	t1, KECCAK1600_NUM_ROUNDS
	blt	a6, t1, 1b

	/* Restore stack */
	ld	s5, 200(sp)
	ld	s4, 208(sp)
	ld	s3, 216(sp)
	ld	s2, 224(sp)
	ld	s1, 232(sp)
	addi	sp, sp, 240
	ret
.endfunc

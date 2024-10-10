// SPDX-License-Identifier: BSD-2-Clause
/*
 * Keccak-f[1600] RV64ID Implementation - State permutation
 * Copyright (C) 2024 Nick Kossifidis <mick@ics.forth.gr>
 */


/*
 * This is the equivalent of keccakf1600_inplaceur
 * using fp registers to hold the state, reducing
 * memory accesses. This is way slower compared to
 * the C implementation but it can be useful in bare
 * metal scenarios where the whole state may be kept
 * in registers, with custom absorb/squeeze functions.
 * In such a case there is no need to save/restore
 * fp registers (unless the rest of the code uses
 * them), nor save/restore state (since it 'll
 * remain in the fp registers for as long as the
 * input message is processed), and the code will
 * also be smaller.
 *
 *
 * a0 -> Pointer to A
 * f0:24 -> Loaded state from A
 * a6 -> Round counter
 * a7 -> Pointer to round_constants
 * s1 - s5 -> C[] / D[]
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

.macro GET_LANE _out, _num
	fmv.x.d	\_out, f\_num
.endm

.macro SET_LANE _num, _in
	fmv.d.x	f\_num, \_in
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

/* Note that a custom absorb/squeeze function
 * could use fp regs directly for accessing
 * the state, in which case we wouldn't need
 * to use any extra memory at all.
 */
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
 * C[i] = A[i] ^ A[i + 5] ^ A[i + 10] ^ A[i + 15] ^ A[i + 20]
 */
.macro COLUMN_PARITY _out, _a, _b, _c, _d, _e
	GET_LANE t0, \_a
	GET_LANE t1, \_b
	GET_LANE t2, \_c
	GET_LANE t3, \_d
	GET_LANE t4, \_e

	xor	t0, t0, t1
	xor	t2, t2, t3
	xor	t4, t4, t0
	xor	\_out, t2, t4
.endm

/**********************\
* KECCAK STEP MAPPINGS *
\**********************/

/*
 * A combined theta_rho step, in-place
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
 */
.macro CHI _a _b _c _d _e
	GET_LANE t1, \_a
	GET_LANE t2, \_b
	GET_LANE t3, \_c
	GET_LANE t4, \_d
	GET_LANE t5, \_e

	_ANDN	a1, t2, t3
	_ANDN   a2, t3, t4
	_ANDN   a3, t4, t5
	_ANDN   a4, t5, t1
	_ANDN   a5, t1, t2

	xor	a1, a1, t1
	xor	a2, a2, t2
	xor	a3, a3, t3
	xor	a4, a4, t4
	xor	a5, a5, t5

	SET_LANE \_a, a1
	SET_LANE \_b, a2
	SET_LANE \_c, a3
	SET_LANE \_d, a4
	SET_LANE \_e, a5
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
	GET_LANE t1, 0
	xor	t0, t0, t1
	SET_LANE 0, t0
.endm

/*************\
* ENTRY POINT *
\*************/

.align 3
.func keccakf1600_state_permute_inplaceur_rv64id
.global keccakf1600_state_permute_inplaceur_rv64id
keccakf1600_state_permute_inplaceur_rv64id:
	/* Save s and fp registers on the stack. No
	 * need to save ra since we won't be calling
	 * any functions from here, and also ignore
	 * the frame pointer. */
	addi	sp, sp, -176
	sd	s1, 168(sp)
	sd	s2, 160(sp)
	sd	s3, 152(sp)
	sd	s4, 144(sp)
	sd	s5, 136(sp)
	SAVE_FP_REGS 0  // fs0:8, fa0:7 -> 0 - 136

	/* Initialize state */
	mv	a6, zero
	la	a7, round_constants

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

	/* Apply theta-rho-pi in-place */

	/* A[0] ^= D[0] (it's unaffected by rho/pi) */
	GET_LANE	t0, 0
	xor		t0, t0, s5
	SET_LANE	0, t0

	GET_LANE	t2, 1	/* Save A[1] for the last step */
	GET_LANE	a2, 6
	GET_LANE	a3, 9
	GET_LANE	a4, 22
	GET_LANE	a5, 14
	THETA_RHO_STEP	a1, a2, s1, 44
	THETA_RHO_STEP	a2, a3, s4, 20
	THETA_RHO_STEP	a3, a4, s2, 61
	THETA_RHO_STEP	a4, a5, s4, 39
	SET_LANE	1, a1
	SET_LANE	6, a2
	SET_LANE	9, a3
	SET_LANE	22, a4

	GET_LANE	a1, 20
	GET_LANE	a2, 2
	GET_LANE	a3, 12
	GET_LANE	a4, 13
	GET_LANE	a5, 19
	THETA_RHO_STEP	t3, a1, s5, 18
	THETA_RHO_STEP	a1, a2, s2, 62
	THETA_RHO_STEP	a2, a3, s2, 43
	THETA_RHO_STEP	a3, a4, s3, 25
	THETA_RHO_STEP	a4, a5, s4, 8
	SET_LANE	14, t3
	SET_LANE	20, a1
	SET_LANE	2, a2
	SET_LANE	12, a3
	SET_LANE	13, a4

	GET_LANE	a1, 23
	GET_LANE	a2, 15
	GET_LANE	a3, 4
	GET_LANE	a4, 24
	GET_LANE	a5, 21
	THETA_RHO_STEP	t3, a1, s3, 56
	THETA_RHO_STEP	a1, a2, s5, 41
	THETA_RHO_STEP	a2, a3, s4, 27
	THETA_RHO_STEP	a3, a4, s4, 14
	THETA_RHO_STEP	a4, a5, s1, 2
	SET_LANE	19, t3
	SET_LANE	23, a1
	SET_LANE	15, a2
	SET_LANE	4, a3
	SET_LANE	24, a4

	GET_LANE	a1, 8
	GET_LANE	a2, 16
	GET_LANE	a3, 5
	GET_LANE	a4, 3
	GET_LANE	a5, 18
	THETA_RHO_STEP	t3, a1, s3, 55
	THETA_RHO_STEP	a1, a2, s1, 45
	THETA_RHO_STEP	a2, a3, s5, 36
	THETA_RHO_STEP	a3, a4, s3, 28
	THETA_RHO_STEP	a4, a5, s3, 21
	SET_LANE	21, t3
	SET_LANE	8, a1
	SET_LANE	16, a2
	SET_LANE	5, a3
	SET_LANE	3, a4

	GET_LANE	a1, 17
	GET_LANE	a2, 11
	GET_LANE	a3, 7
	GET_LANE	a4, 10
	THETA_RHO_STEP	t3, a1, s2, 15
	THETA_RHO_STEP	a1, a2, s1, 10
	THETA_RHO_STEP	a2, a3, s2, 6
	THETA_RHO_STEP	a3, a4, s5, 3
	THETA_RHO_STEP	a4, t2, s1, 1
	SET_LANE	18, t3
	SET_LANE	17, a1
	SET_LANE	11, a2
	SET_LANE	7, a3
	SET_LANE	10, a4

	CHI 0 1 2 3 4
	CHI 5 6 7 8 9
	CHI 10 11 12 13 14
	CHI 15 16 17 18 19
	CHI 20 21 22 23 24

	IOTA

	addi	a6, a6, 1
	li	t1, KECCAK1600_NUM_ROUNDS
	blt	a6, t1, 1b
	STORE_STATE

	RESTORE_FP_REGS 0
	ld	s5, 136(sp)
	ld	s4, 144(sp)
	ld	s3, 152(sp)
	ld	s2, 160(sp)
	ld	s1, 168(sp)
	addi	sp, sp, 176
	ret
.endfunc

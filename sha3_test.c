// SPDX-License-Identifier: BSD-2-Clause
/*
 * SHA3 C / RV64 Implementation - Test suite
 * Copyright (C) 2019 Nick Kossifidis <mick@ics.forth.gr>
 */

#include <stdio.h>	/* For printf() */
#include <stdlib.h>	/* For malloc() */
#include <time.h>	/* For clock() */
#include <string.h>	/* For memset() */
#include "keccak1600.h"
#include "sha3.h"

/*
 * To verify the output check out:
 * https://www.di-mgt.com.au/sha_testvectors.html
 */

static void
sha3_print(const char* md, size_t md_len)
{
	int i = 0;
	for(i = 0; i < md_len; i++)
		printf("%02X", md[i] & 0xFF);
	printf("\n");
}

static clock_t
sha3_test(int print, char* amillion_as) {
	char md256[32] = {0};
	char md512[64] = {0};
	clock_t start = 0;
	clock_t end = 0;

	start = clock();

	sha3_256_oneshot("", 0, md256);
	if(print) {
		printf("SHA3-256 of empty string:\t");
		sha3_print((const char*) md256, 32);
	}

	sha3_512_oneshot("", 0, md512);
	if(print) {
		printf("SHA3-512 of empty string:\t");
		sha3_print((const char*) md512, 64);
	}

	sha3_256_oneshot("abc", 3, md256);
	if(print) {
		printf("SHA3-256 of \"abc\":\t\t");
		sha3_print((const char*) md256, 32);
	}

	sha3_512_oneshot("abc", 3, md512);
	if(print) {
		printf("SHA3-512 of \"abc\":\t\t");
		sha3_print((const char*) md512, 64);
	}

	sha3_256_oneshot("test", 4, md256);
	if(print) {
		printf("SHA3-256 of \"test\":\t\t");
		sha3_print((const char*) md256, 32);
	}

	sha3_512_oneshot("test", 3, md512);
	if(print) {
		printf("SHA3-512 of \"test\":\t\t");
		sha3_print((const char*) md512, 64);
	}

	sha3_256_oneshot(amillion_as, 1000000, md256);
	if(print) {
		printf("SHA3-256 of 1mil 'a's:\t\t");
		sha3_print((const char*) md256, 32);
	}

	sha3_512_oneshot(amillion_as, 1000000, md512);
	if(print) {
		printf("SHA3-512 of 1mil 'a's:\t\t");
		sha3_print((const char*) md512, 64);
	}

	end = clock();

	return (end - start);
}

int
main()
{
	double test_dur = 0;
	double ema = 0;
	double n = 10;
	char *amillion_as = NULL;
	int i = 0;

	amillion_as = malloc(sizeof(char) * 1000000);
	memset(amillion_as, 0x61, sizeof(char) * 1000000);

#ifdef OSSL_BUILD
	printf("\nOpenSSL implementation\n");
	printf("======================\n");
	for(i = 0; i < n; i++) {
		test_dur = (double) sha3_test(!i, amillion_as);
		ema = (test_dur + (n - 1) * ema) / n;
	}
	printf("Test took an avg of %lg sec (%lg clock ticks)\n", ema / CLOCKS_PER_SEC, ema);
#else
	printf("\nSimple implementation\n");
	printf("=====================\n");
	keccakf1600_set_permutation_function(&keccakf1600_state_permute_simple, 0);
	for(i = 0; i < n; i++) {
		test_dur = (double) sha3_test(!i, amillion_as);
		ema = (test_dur + (n - 1) * ema) / n;
	}
	printf("Test took an avg of %lg sec (%lg clock ticks)\n", ema / CLOCKS_PER_SEC, ema);
	ema = 0;

	printf("\nCompact implementation\n");
	printf("======================\n");
	keccakf1600_set_permutation_function(&keccakf1600_state_permute_compact, 0);
	for(i = 0; i < n; i++) {
		test_dur = (double) sha3_test(!i, amillion_as);
		ema = (test_dur + (n - 1) * ema) / n;
	}
	printf("Test took an avg of %lg sec (%lg clock ticks)\n", ema / CLOCKS_PER_SEC, ema);
	ema = 0;

	printf("\nIn-place unrolled\n");
	printf("=================\n");
	keccakf1600_set_permutation_function(&keccakf1600_state_permute_inplaceur, 0);
	for(i = 0; i < n; i++) {
		test_dur = (double) sha3_test(!i, amillion_as);
		ema = (test_dur + (n - 1) * ema) / n;
	}
	printf("Test took an avg of %lg sec (%lg clock ticks)\n", ema / CLOCKS_PER_SEC, ema);
	ema = 0;

	printf("\nUnrolled with intermediate state (cache friendly)\n");
	printf("=================================================\n");
	keccakf1600_set_permutation_function(&keccakf1600_state_permute_intermediateur, 0);
	for(i = 0; i < n; i++) {
		test_dur = (double) sha3_test(!i, amillion_as);
		ema = (test_dur + (n - 1) * ema) / n;
	}
	printf("Test took an avg of %lg sec (%lg clock ticks)\n", ema / CLOCKS_PER_SEC, ema);
	ema = 0;

	printf("\nUnrolled with intermediate state + early parity\n");
	printf("===============================================\n");
	keccakf1600_set_permutation_function(&keccakf1600_state_permute_intermediateur_ep, 0);
	for(i = 0; i < n; i++) {
		test_dur = (double) sha3_test(!i, amillion_as);
		ema = (test_dur + (n - 1) * ema) / n;
	}
	printf("Test took an avg of %lg sec (%lg clock ticks)\n", ema / CLOCKS_PER_SEC, ema);
	ema = 0;

	printf("\nUnrolled with intermediate state + lane complementing\n");
	printf("=====================================================\n");
	keccakf1600_set_permutation_function(&keccakf1600_state_permute_intermediateur_lc, 1);
	for(i = 0; i < n; i++) {
		test_dur = (double) sha3_test(!i, amillion_as);
		ema = (test_dur + (n - 1) * ema) / n;
	}
	printf("Test took an avg of %lg sec (%lg clock ticks)\n", ema / CLOCKS_PER_SEC, ema);
	ema = 0;
#ifdef RVASM_IMPL
	printf("\nUnrolled with intermediate state (RV64I)\n");
	printf("========================================\n");
	keccakf1600_set_permutation_function(&keccakf1600_state_permute_rv64i, 0);
	for(i = 0; i < n; i++) {
		test_dur = (double) sha3_test(!i, amillion_as);
		ema = (test_dur + (n - 1) * ema) / n;
	}
	printf("Test took an avg of %lg sec (%lg clock ticks)\n", ema / CLOCKS_PER_SEC, ema);
	ema = 0;

	printf("\nIn-place unrolled (RV64ID)\n");
	printf("==========================\n");
	keccakf1600_set_permutation_function(&keccakf1600_state_permute_rv64id, 0);
	for(i = 0; i < n; i++) {
		test_dur = (double) sha3_test(!i, amillion_as);
		ema = (test_dur + (n - 1) * ema) / n;
	}
	printf("Test took an avg of %lg sec (%lg clock ticks)\n", ema / CLOCKS_PER_SEC, ema);
#endif /* RVASM_IMPL */
#endif /* OSSL_BUILD */
}

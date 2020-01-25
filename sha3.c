// SPDX-License-Identifier: BSD-2-Clause
/*
 * SHA3 C / RV64 Implementation
 * Copyright (C) 2019 Nick Kossifidis <mick@ics.forth.gr>
 */

#include "keccak1600.h"
#include "sha3.h"

static void
sha3_oneshot(const void *msg, size_t msg_len, void *md, size_t md_len)
{
	keccakf1600_oneshot(msg, msg_len, md, md_len, 0x06);
}

/**************\
* ENTRY POINTS *
\**************/

void sha3_256_oneshot(const void *msg, size_t msg_len, void *md)
{
	sha3_oneshot(msg, msg_len, md, 32);
}

void sha3_512_oneshot(const void *msg, size_t msg_len, void *md)
{
	sha3_oneshot(msg, msg_len, md, 64);
}

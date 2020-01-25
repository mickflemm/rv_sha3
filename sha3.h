// SPDX-License-Identifier: BSD-2-Clause
/*
 * SHA3 C / RV64 Implementation
 * Copyright (C) 2019 Nick Kossifidis <mick@ics.forth.gr>
 */

#include <stddef.h>	/* For size_t */

void sha3_256_oneshot(const void *msg, size_t msg_len, void* md);
void sha3_512_oneshot(const void *msg, size_t msg_len, void* md);

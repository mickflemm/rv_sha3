// SPDX-License-Identifier: BSD-2-Clause
/*
 * SHA3 OpenSSL Wrapper for testing/comparison
 * Copyright (C) 2019 Nick Kossifidis <mick@ics.forth.gr>
 */

#include <openssl/evp.h>
#include "sha3.h"


/**************\
* ENTRY POINTS *
\**************/

void
sha3_256_oneshot(const void *msg, size_t msg_len, void* md)
{
	EVP_MD_CTX *ctx = EVP_MD_CTX_create();
	unsigned int md_len = 32;
	EVP_DigestInit_ex(ctx, EVP_sha3_256(), NULL);
	EVP_DigestUpdate(ctx, msg, msg_len);
	EVP_DigestFinal_ex(ctx, md, &md_len);
	EVP_MD_CTX_destroy(ctx);
}

void
sha3_512_oneshot(const void *msg, size_t msg_len, void* md)
{
	EVP_MD_CTX *ctx = EVP_MD_CTX_create();
	unsigned int md_len = 64;
	EVP_DigestInit_ex(ctx, EVP_sha3_512(), NULL);
	EVP_DigestUpdate(ctx, msg, msg_len);
	EVP_DigestFinal_ex(ctx, md, &md_len);
	EVP_MD_CTX_destroy(ctx);
}

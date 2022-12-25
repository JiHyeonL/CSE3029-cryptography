/*
 * Copyright 2020-2022. Heekuck Oh, all rights reserved
 * 이 프로그램은 한양대학교 ERICA 소프트웨어학부 재학생을 위한 교육용으로 제작되었다.
 */
#ifdef __linux__
#include <bsd/stdlib.h>
#elif __APPLE__
#include <stdlib.h>
#else
#include <stdlib.h>
#endif
#include <string.h>
#include <gmp.h>
#include "pkcs.h"
#include "sha2.h"

/*
 * rsa_generate_key() - generates RSA keys e, d and n in octet strings.
 * If mode = 0, then e = 65537 is used. Otherwise e will be randomly selected.
 * Carmichael's totient function Lambda(n) is used.
 */
void rsa_generate_key(void *_e, void *_d, void *_n, int mode)
{
    mpz_t p, q, lambda, e, d, n, gcd;
    gmp_randstate_t state;
    
    /*
     * Initialize mpz variables
     */
    mpz_inits(p, q, lambda, e, d, n, gcd, NULL);
    gmp_randinit_default(state);
    gmp_randseed_ui(state, arc4random());
    /*
     * Generate prime p and q such that 2^(RSAKEYSIZE-1) <= p*q < 2^RSAKEYSIZE
     */
    do {
        do {
            mpz_urandomb(p, state, RSAKEYSIZE/2);
            mpz_setbit(p, 0);
            mpz_setbit(p, RSAKEYSIZE/2-1);
        } while (mpz_probab_prime_p(p, 50) == 0);
        do {
            mpz_urandomb(q, state, RSAKEYSIZE/2);
            mpz_setbit(q, 0);
            mpz_setbit(q, RSAKEYSIZE/2-1);
        } while (mpz_probab_prime_p(q, 50) == 0);
        mpz_mul(n, p, q);
    } while (!mpz_tstbit(n, RSAKEYSIZE-1));
    /*
     * Generate e and d using Lambda(n)
     */
    mpz_sub_ui(p, p, 1);
    mpz_sub_ui(q, q, 1);
    mpz_lcm(lambda, p, q);
    if (mode == 0)
        mpz_set_ui(e, 65537);
    else do {
        mpz_urandomb(e, state, RSAKEYSIZE);
        mpz_gcd(gcd, e, lambda);
    } while (mpz_cmp(e, lambda) >= 0 || mpz_cmp_ui(gcd, 1) != 0);
    mpz_invert(d, e, lambda);
    /*
     * Convert mpz_t values into octet strings
     */
    mpz_export(_e, NULL, 1, RSAKEYSIZE/8, 1, 0, e);
    mpz_export(_d, NULL, 1, RSAKEYSIZE/8, 1, 0, d);
    mpz_export(_n, NULL, 1, RSAKEYSIZE/8, 1, 0, n);
    /*
     * Free the space occupied by mpz variables
     */
    mpz_clears(p, q, lambda, e, d, n, gcd, NULL);
}

/*
 * rsa_cipher() - compute m^k mod n
 * If m >= n then returns PKCS_MSG_OUT_OF_RANGE, otherwise returns 0 for success.
 */
static int rsa_cipher(void *_m, const void *_k, const void *_n)
{
    mpz_t m, k, n;
    
    /*
     * Initialize mpz variables
     */
    mpz_inits(m, k, n, NULL);
    /*
     * Convert big-endian octets into mpz_t values
     */
    mpz_import(m, RSAKEYSIZE/8, 1, 1, 1, 0, _m);
    mpz_import(k, RSAKEYSIZE/8, 1, 1, 1, 0, _k);
    mpz_import(n, RSAKEYSIZE/8, 1, 1, 1, 0, _n);
    /*
     * Compute m^k mod n
     */
    if (mpz_cmp(m, n) >= 0) {
        mpz_clears(m, k, n, NULL);
        return PKCS_MSG_OUT_OF_RANGE;
    }
    mpz_powm(m, m, k, n);
    /*
     * Convert mpz_t m into the octet string _m
     */
    mpz_export(_m, NULL, 1, RSAKEYSIZE/8, 1, 0, m);
    /*
     * Free the space occupied by mpz variables
     */
    mpz_clears(m, k, n, NULL);
    return 0;
}

/* sel_shaFunc() - 사용할 sha2 함수 선택해서 호출 */
void sel_shaFunc(const unsigned char *message, unsigned int len, unsigned char *digest, int sha2_ndx)
{
	switch (sha2_ndx)
	{
		case SHA224:
			sha224(message, len, digest);
			break;
		case SHA256:
			sha256(message, len, digest);
			break;
		case SHA384:
			sha384(message, len, digest);
			break;
		case SHA512:
			sha512(message, len, digest);
			break;
		case SHA512_224:
			sha512_224(message, len, digest);
			break;
		case SHA512_256:
			sha512_256(message, len, digest);
			break;
	}
}

/* sel_shaLen() - 사용할 sha2 함수 output 길이 return */
int sel_shaLen(int sha2_ndx)
{
	switch (sha2_ndx)
	{
		case SHA224:
			return SHA224_DIGEST_SIZE;
		case SHA256:
			return SHA256_DIGEST_SIZE;
		case SHA384:
			return SHA384_DIGEST_SIZE;
		case SHA512:
			return SHA512_DIGEST_SIZE;
		case SHA512_224:
			return SHA224_DIGEST_SIZE;
		case SHA512_256:
			return SHA256_DIGEST_SIZE;
	}
	return 0;
}

/* i2osp() - int x를 octet string으로 바꿔준다.(mgf에서 사용) */
void i2osp(int x, int xLen, unsigned char *digit)
{
	for (int j = 0; j < xLen; j++)
	{
		digit[3 - j] = x & 0x000000ff;
		x >>= 8;
	}
}

/* mgf() - mgfSeed를 maskLen 길이인 mask로 return */
unsigned char *mgf(const unsigned char *mgfSeed, size_t seedLen, unsigned char *mask, size_t maskLen, int sha2_ndx)
{
	uint32_t counter = maskLen / seedLen + (maskLen % seedLen ? 1 : 0);
	unsigned char mgfStore[seedLen + 4], t[counter * seedLen];

	memcpy(mgfStore, mgfSeed, seedLen);	// mgfSeed를 seedLen만큼 mgfStore에 copy

	for (int i = 0; i < counter; i++) {
		i2osp(i, 4, mgfStore+seedLen);	// counter를 길이가 4인 octet string으로 변환 후 mgfStore의 seedLen 뒤에 저장
		sel_shaFunc(mgfStore, seedLen + 4, t + i * seedLen, sha2_ndx);	// mgfStore를 hash한 값을 t의 i * seedLen 뒤에 저장
	}

	memcpy(mask, t, maskLen);	// t를 길이가 maskLen인 mask로 copy

	return mask;
}

/*
 * rsaes_oaep_encrypt() - RSA encrytion with the EME-OAEP encoding method
 * 길이가 len 바이트인 메시지 m을 공개키 (e,n)으로 암호화한 결과를 c에 저장한다.
 * label은 데이터를 식별하기 위한 라벨 문자열로 NULL을 입력하여 생략할 수 있다.
 * sha2_ndx는 사용할 SHA-2 해시함수 색인 값으로 SHA224, SHA256, SHA384, SHA512,
 * SHA512_224, SHA512_256 중에서 선택한다. c의 크기는 RSAKEYSIZE와 같아야 한다.
 * 성공하면 0, 그렇지 않으면 오류 코드를 넘겨준다.
 */
int rsaes_oaep_encrypt(const void *m, size_t mLen, const void *label, const void *e, const void *n, void *c, int sha2_ndx)
{
    /* 목표: EM을 공개키 (e,n)을 사용해 EM^e mod n 으로 암호화 */

    int hlen = sel_shaLen(sha2_ndx);
    int keylen = RSAKEYSIZE/8;
    int dblen = keylen-hlen-1;
    int pslen;
    unsigned char hlabel[hlen]; 
    unsigned char db[dblen], seed[hlen], dbMask[dblen], maskedDB[dblen], seedMask[hlen], maskedSeed[hlen];
    unsigned char em[keylen];

    /* 1. length checking */
    if (mLen > keylen - 2*hlen -2)
        return 2;   // 2) PKCS_MSG_TOO_LONG

    /* 2. EME-OAEP encoding */
    sel_shaFunc(label, strlen(label), hlabel, sha2_ndx);
    pslen = keylen-mLen-2*hlen-2;
    memcpy(db, hlabel, hlen);
    memset(db+hlen, 0x00, pslen);
    memset(db+hlen+pslen, 0x01, 1);
    memcpy(db+hlen+pslen+1, m, mLen);   // db 완성(길이: RSAKEYSIZE/8-hlen-1)
    

    arc4random_buf(seed,hlen); 
    mgf(seed, hlen, dbMask, dblen, sha2_ndx);    // mgf1

    for(int i = 0; i < dblen; i++) {
        maskedDB[i] = db[i] ^ dbMask[i];    // maskedDB = db xor mgf1
    }

    
    mgf(maskedDB, dblen, seedMask, hlen, sha2_ndx);	// mgf2
    
    for(int i = 0; i < hlen; i++) {
        maskedSeed[i] = seed[i] ^ seedMask[i];  // maskedSeed = seed xor mgf2
    }

    memset(em, 0x00, 1);
    memcpy(em+1, maskedSeed, hlen);
    memcpy(em+1+hlen, maskedDB, dblen);  // em 완성(길이: RSAKEYSIZE/8)


    /* 3. RSA encryption */
    if(rsa_cipher((void *)em, e, n) == 1)
        return 1;   // 1) PKCS_MSG_OUT_OF_RANGE 

    memcpy(c, (void *)em, keylen);

    return 0;
}

/*
 * rsaes_oaep_decrypt() - RSA decrytion with the EME-OAEP encoding method
 * 암호문 c를 개인키 (d,n)을 사용하여 원본 메시지 m과 길이 len을 회복한다.
 * label과 sha2_ndx는 암호화할 때 사용한 것과 일치해야 한다.
 * 성공하면 0, 그렇지 않으면 오류 코드를 넘겨준다.
 */

int rsaes_oaep_decrypt(void *m, size_t *mLen, const void *label, const void *d, const void *n, const void *c, int sha2_ndx)
{
    int hLen = sel_shaLen(sha2_ndx);
    int DBLen = (RSAKEYSIZE/8) - 1 - hLen;
    int i = 0, PSLen = 0;
    unsigned char *_label, *EM;
    unsigned char seed[hLen], DB[DBLen], lHash[hLen], lHashPrime[hLen], mask[DBLen];

    /* 암호문 c RSA 복호화 = EM */
    rsa_cipher((void *)c, d, n);
    EM = (unsigned char *) c;

    /* sha2_ndx에서 정한 해시 함수로 label의 해시값 구하고 lHash에 저장 */
    _label = (unsigned char *) label;
    sel_shaFunc(_label, strlen((char *)label), lHash, sha2_ndx);

    /* EM = 0x00||maskedSeed||maskedDB */
    if (EM[0])  // EM의 첫번째 옥텟이 0x00인지 확인
        return PKCS_INITIAL_NONZERO;
    memcpy(seed, EM + 1, hLen);
    memcpy(DB, EM + 1 + hLen, DBLen);


    /* seedMask = MGF(maskedDB, hLen) */
    mgf(DB, DBLen, mask, hLen, sha2_ndx);


    /* seed = maskedSeed XOR seedMask */
    for (i = 0; i < hLen; i++)
        seed[i] ^= mask[i];

    /* DBMask = MGF(seed, DBLen) */
    mgf(seed, hLen, mask, DBLen, sha2_ndx);

    /* DB = maskedDB XOR DBMask */
    for (i = 0; i < DBLen; i++)
        DB[i] ^= mask[i];
   
    /* DB = lHashPrime||PS||0x01||M */
    memcpy(lHashPrime, DB, hLen);

    /* lHashPrime = lHash인지 확인 */
    for (i = 0; i < hLen; i++)
        if (lHash[i] != lHashPrime[i])
            return PKCS_HASH_MISMATCH;

    /* PS 길이 세기 */
    for (i = hLen; i < DBLen; i++)
    {
        if (DB[i] != 0x00)
            break;
        PSLen++;
    }

    /* PS와 원래 메세지를 구분하는 옥텟 "0x01"이 있는지 확인  */
    if (DB[hLen + PSLen] != 0x01)
        return PKCS_INVALID_PS;

    /* 원래 메세지의 길이 계산 및 저장 */
    *mLen = DBLen - hLen - PSLen - 1;

    /* m에 원래 메세지 저장 */
    memcpy(m, DB + hLen + PSLen + 1, (*mLen));

    return 0;
}

/*
 * rsassa_pss_sign - RSA Signature Scheme with Appendix
 * 길이가 len 바이트인 메시지 m을 개인키 (d,n)으로 서명한 결과를 s에 저장한다.
 * s의 크기는 RSAKEYSIZE와 같아야 한다. 성공하면 0, 그렇지 않으면 오류 코드를 넘겨준다.
 */
int rsassa_pss_sign(const void *m, size_t mLen, const void *d, const void *n, void *s, int sha2_ndx)
{
    /* 메시지 길이 검사 */
    if (mLen > 0x1fffffffffffffff) return PKCS_MSG_TOO_LONG;
    
    /* 해시함수에 따른 사이즈 설정 (바이트 단위) */ 
    int mSize = sel_shaLen(sha2_ndx); // 메시지 길이
    int KEYSIZE = RSAKEYSIZE / 8; // 키 길이
    int dbSize = KEYSIZE - (mSize+1); // DB 길이
    int psSize = KEYSIZE-(2*(mSize+1)); // PS 길이

    /* 메모리 할당 및 변수 생성 */ 
    unsigned char mHash[mSize];
    unsigned char mPrime[8+(mSize*2)];
    unsigned char H[mSize];
    unsigned char DB[dbSize];
    unsigned char MGF[dbSize];
    unsigned char MDB[dbSize];
    unsigned char EM[KEYSIZE];
    uint8_t salt[mSize];
    int res; // 결과값
    
    /* mHash 생성 */
    sel_shaFunc((unsigned char*)m, mLen, mHash, sha2_ndx);
    
    /* salt 생성 */
    arc4random_buf(&salt, sizeof(uint8_t)*((mSize)));

    /* 결합 (M') */    
    memset(mPrime, 0x00, sizeof(char)*8); // M' = M' || header(8bytes 0x00)
    memcpy(mPrime+sizeof(char)*8, mHash, sizeof(char)*mSize); // M' = M' || mHash
    memcpy(mPrime+sizeof(char)*(8+mSize), salt, sizeof(char)*mSize); // M' = M' || salt

    /* H 생성 */
    sel_shaFunc(mPrime, sizeof(mPrime), H, sha2_ndx);

    /* 결합 (DB) */
    memset(DB, 0x00, sizeof(char)*psSize); // DB = PS
    memset(DB+sizeof(char)*psSize, 0x01, sizeof(char)); // DB = DB || '0x01'
    memcpy(DB+sizeof(char)*(psSize+1), salt, sizeof(char)*mSize); // DB = DB || salt

    /* MGF 계산 */
    mgf(H, mSize, MGF, sizeof(MGF), sha2_ndx);

    /* maskedDB 생성 */
    for (int i = 0; i < dbSize; ++i) { // DB xor MGF
        MDB[i] = DB[i] ^ MGF[i];
    }
    if(MDB[0] >> 7) MDB[0] &= 0x7f; // MDB의 가장 왼쪽값이 1일경우 0으로 변경

    /* 결합 (EM) */
    memcpy(EM, MDB, sizeof(char)*dbSize); // EM = MDB
    memcpy(EM+sizeof(char)*dbSize, H, sizeof(char)*mSize); // EM = EM || H
    memset(EM+sizeof(char)*(dbSize+mSize), 0xBC, sizeof(char)); // EM = EM || '0xBC'

    /* 암호화 */
    res = rsa_cipher(EM, d, n); // 암호화
    if (!res) {
        memcpy(s, EM, sizeof(EM)); // 암호화 결과 저장
    }
    return res; // 결과값 반환
}

/*
 * rsassa_pss_verify - RSA Signature Scheme with Appendix
 * 길이가 len 바이트인 메시지 m에 대한 서명이 s가 맞는지 공개키 (e,n)으로 검증한다.
 * 성공하면 0, 그렇지 않으면 오류 코드를 넘겨준다.
 */
int rsassa_pss_verify(const void *m, size_t mLen, const void *e, const void *n, const void *s, int sha2_ndx)
{
    int sha_size=sel_shaLen(sha2_ndx);
    int key_size=RSAKEYSIZE/8;
    int db_size=key_size-sha_size-1;
    int mpr_size=sha_size*2+8;
    unsigned char em[key_size];
    unsigned char mHash[sha_size];
    unsigned char maskedDb[db_size];
    unsigned char h[sha_size];
    unsigned char dbMask[db_size];
    unsigned char db[db_size];
    unsigned char dbPad[db_size-sha_size];
    unsigned char salt[sha_size];
    unsigned char mPrime[mpr_size];
    unsigned char hashPrime[sha_size];

    /* m의 길이가 hash function의 제한(2^61-1 octets for SHA-1)보다 클 경우 - PKCS_MSG_TOO_LONG */
    if (mLen>0x1fffffffffffffff) return PKCS_MSG_TOO_LONG;
    
    /* mHash 만들기 */
    sel_shaFunc(m, mLen, mHash, sha2_ndx);

    /* em 구하기 */
    memcpy(em, s, key_size);
    rsa_cipher((void *)em, e, n);

    /* em의 마지막 octet이 0xbc가 아닐 경우 - PKCS_INVALID_LAST */
    if(em[key_size-1]!=0xbc) return PKCS_INVALID_LAST;

    /* maskDb(leftmost emLen - hLen - 1 octets of EM) 구하기 */
    memcpy(maskedDb, em, db_size);

    /* h(next hLen octets) 구하기 */
    memcpy(h, em+db_size, sha_size);

    /* em의 첫 번째 bit가 0이 아닐 경우 - PKCS_INVALID_INIT */
    if((em[0]>>7)!=0) return PKCS_INVALID_INIT;

    /* dbMask = MGF(H, emLen-hLen-1) 구하기 */
    mgf(h, sha_size, dbMask, db_size, sha2_ndx);

    /* db = maskedDb xor dbMask */
    for (int i=0; i<db_size; i++) {
        db[i]=maskedDb[i]^dbMask[i];
    }

    /* db padding 구하기 */
    memcpy(dbPad, db, db_size-sha_size);

    /* db padding이 0000....0001 꼴이 아닐 경우 - PKCS_INVALID_PD2 */
    if(dbPad[db_size-sha_size-1]!=0x01) return PKCS_INVALID_PD2;
    
    /* salt 구하기  */
    memcpy(salt, db+db_size-sha_size, sha_size);

    /* m' 구하기 ((0x)00 00 00 00 00 00 00 00 || mHash || salt) */
    memset(mPrime, 0x0, 8);
    memcpy(mPrime+8, mHash, sha_size);
    memcpy(mPrime+8+sha_size, salt, sha_size);

    /* H' = hash(m') */
    sel_shaFunc(mPrime, mpr_size, hashPrime, sha2_ndx);

    /* if(H!=H') return PKCS_HASH_MISMATCH */
    if(memcmp(hashPrime, h, sha_size)!=0) return PKCS_HASH_MISMATCH;
    return 0;
}

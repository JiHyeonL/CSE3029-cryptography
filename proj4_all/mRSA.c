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
#include "mRSA.h"

/*
 * mod_add() - computes a + b mod m
 */
static uint64_t mod_add(uint64_t a, uint64_t b, uint64_t m)
{
    a = a % m;
    b = b % m;
    return a >= m-b ? (a-(m-b)) % m : (a+b) % m;
}

/*
 * mod_sub() - computes a - b mod m
 */
static uint64_t mod_sub(uint64_t a, uint64_t b, uint64_t m)
{
    a = a % m;
    b = b % m;
    return a < b ? (a-b+m) % m :(a-b) % m;
}

/*
 * mod_mul() - computes a * b mod m
 */
static uint64_t mod_mul(uint64_t a, uint64_t b, uint64_t m)
{
    uint64_t r = 0;
    while (b > 0) {
        if (b & 1)
            r = mod_add(r, a, m);
        b = b >> 1;
        a = mod_add(a, a, m);
    }
    return r;
}

/*
 * mod_pow() - computes a^b mod m
 */
static uint64_t mod_pow(uint64_t a, uint64_t b, uint64_t m)
{
    uint64_t r = 1;
    while (b > 0) {
        if (b & 1)
            r = mod_mul(r, a, m);
        b = b >> 1;
        a = mod_mul(a, a, m);
    }
    return r;    
}

/*
 * gcd() - Euclidean algorithm
 */
static uint64_t gcd(uint64_t a, uint64_t b)
{
    uint64_t temp = 0;
    while(1) {
        if(a == 0 || b == 0) {
            return a + b; //a나 b가 0이면 a나 b가 최대공약수
        }

        temp = a;
        a = b; // a = b
        b = temp % b; // b = a % b
    }    
}

/*
 * mul_inv() - computes multiplicative inverse a^-1 mod m
 * It returns 0 if no inverse exist.
 */
static uint64_t mul_inv(uint64_t a, uint64_t m)
{
    uint64_t d0 = a, d1 = m;
    uint64_t x0 = 1, x1 = 0, q, temp;

    while(d1 > 1){
        q = d0 / d1;

        temp = mod_sub(d0 ,mod_mul(q, d1, m), m); // d0 – q * d1
        d0 = d1;
        d1 = temp;

        temp = mod_sub(x0, mod_mul(q, x1, m), m); // x0 – q * x1
        x0 = x1;
        x1 = temp;
    }

    if(d1 == 1) // 역이 존재할 경우
        return (x1 > 0 ? x1 : x1+m); // x1이 음수일 경우 x1+m return
    else
        return 0;
}

/*
 * Miller-Rabin Primality Testing against small sets of bases
 *
 * if n < 2^64,
 * it is enough to test a = 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, and 37.
 *
 * if n < 3317044064679887385961981,
 * it is enough to test a = 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, and 41.
 */
static const uint64_t a[BASELEN] = {2,3,5,7,11,13,17,19,23,29,31,37};

/*
 * miller_rabin() - Miller-Rabin Primality Test (deterministic version)
 *
 * n > 3, an odd integer to be tested for primality
 * It returns 1 if n is prime, 0 otherwise.
 */
static int find_prime(uint64_t n,uint64_t k,uint64_t temp)  // (a[i]^p)^(2^j) mod n 계산
{
    for(int j = 0; j < k; j++) {
        if (temp == n-1) {  // prime
            return 1;
        }
        temp = mod_mul(temp,temp,n); // 이전 temp값을 제곱하면 다음 temp값이 나온다.(a[i]^p, a[i]^2p, a[i]^4p, ...)
    }
    return 0; // composite
}

static int miller_rabin(uint64_t n)
{
    uint64_t k,q,temp;

    if(n == 0 || n == 1 || (n > 2 && n % 2 == 0))
        return COMPOSITE;

    // k,q 찾기
    q = (n - 1);
    k = 0;
    while(q % 2 == 0){
        q = q / 2;
        k = k + 1;
    }

    //PRIME or COMPOSITE
    for(int i = 0; i < BASELEN; i++) {
        if(a[i] >= n-1)
            continue;

        temp = mod_pow(a[i],q,n);
        if (temp == 1)   // prime
            continue;  

        if (find_prime(n,k,temp) == 1)  // prime
            continue;

        return COMPOSITE; // composite
    }
    return PRIME; // prime    
}

/*
 * mRSA_generate_key() - generates mini RSA keys e, d and n
 *
 * Carmichael's totient function Lambda(n) is used.
 */
void mRSA_generate_key(uint64_t *e, uint64_t *d, uint64_t *n)
{
    uint64_t p = 0, q = 0, lamb;

    // 1. p,q 찾기
    while (p * q < MINIMUM_N) {
        while (1) {
            arc4random_buf(&p, sizeof(uint32_t)); // random p값 생성
            if (miller_rabin(p)) // p가 소수이면
                break;
        }
        while (1) {
            arc4random_buf(&q, sizeof(uint32_t)); // random q값 생성
            if (miller_rabin(q)) // q가 소수이면
                break;
        }
    }


    // 2. n = p * q, lamb = (p-1)(q-1) 계산
    *n = p * q; 
    lamb = (p-1)*(q-1)/gcd(p-1,q-1);


    // 3. gcd(e,lamb) = 1이고 1 < e < lamb인 e 찾기
    while (1) {
        arc4random_buf(e,sizeof(uint64_t)); // random e값 생성
        if ((1 < *e) && (*e < lamb) && (gcd(*e,lamb) == 1))
            break;       
    }


    // 4. ed ≡ 1 mod lamb 인 d 찾기
    *d = mul_inv(*e,lamb);    
}

/*
 * mRSA_cipher() - compute m^k mod n
 *
 * If data >= n then returns 1 (error), otherwise 0 (success).
 */
int mRSA_cipher(uint64_t *m, uint64_t k, uint64_t n)
{
    if(*m >= n)
        return 1;
    else {
        *m = mod_pow(*m, k, n);
        return 0;
    }    
}

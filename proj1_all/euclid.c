/*
 * Copyright 2020-2022. Heekuck Oh, all rights reserved
 * 이 프로그램은 한양대학교 ERICA 소프트웨어학부 재학생을 위한 교육용으로 제작되었다.
 */
#include "euclid.h"

/*
 * gcd() - Euclidean algorithm
 *
 * 유클리드 알고리즘 gcd(a,b) = gcd(b,a mod b)를 사용하여 최대공약수를 계산한다.
 * 만일 a가 0이면 b가 최대공약수가 된다. 그 반대도 마찬가지이다.
 * a, b가 모두 음이 아닌 정수라고 가정한다.
 * 재귀함수 호출을 사용하지 말고 while 루프를 사용하여 구현하는 것이 빠르고 좋다.
 */

int gcd(int a, int b)
{
    int temp = 0;
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
 * xgcd() - Extended Euclidean algorithm
 *
 * 확장유클리드 알고리즘은 두 수의 최대공약수 gcd(a,b) = ax + by 식을
 * 만족하는 x와 y를 계산하는 알고리즘이다. 강의노트를 참조하여 구현한다.
 * a, b가 모두 음이 아닌 정수라고 가정한다.   
 */
int xgcd(int a, int b, int *x, int *y)
{
    int d0, d1, x0, x1, y0, y1, q, temp;
    d0 = a, d1 = b, x0 = 1, x1 = 0, y0 = 0, y1 = 1;

    while(1) {
        q = d0 / d1;  // q = d0 div d1
        
        temp = d1;
        d1 = d0 - q*d1; // d2 = d0 – q * d1
        d0 = temp;   

        temp = x1;
        x1 = x0 - q*x1; // x2 = x0 – q * x1
        x0 = temp;

        temp = y1;
        y1 = y0 - q*y1; // y2 = y0 – q * y1
        y0 = temp;

        if(d1 == 0) { // d1이 0이면
            *x = x0; // x값
            *y = y0; // y값
            return d0; // 최대공약수 d0 return
        }
    }
}

/*
 * mul_inv() - computes multiplicative inverse a^-1 mod m
 *
 * 모듈로 m에서 a의 곱의 역인 a^-1 mod m을 구한다.
 * 만일 역이 존재하지 않으면 0을 리턴한다.
 * 확장유클리드 알고리즘을 변형하여 구현한다. 강의노트를 참조한다.
 */
int mul_inv(int a, int m)
{
    int d0 = a, d1 = m;
    int x0 = 1, x1 = 0, q, temp;

    while(d1 > 1) {
        q = d0 / d1; // q = d0 div d1

        d0 = d0 - q * d1; // d2 = d0 – q * d1
        temp = d0;
        d0 = d1;
        d1 = temp;

        x0 = x0 - q*x1; // x2 = x0 – q * x1
        temp = x0;
        x0 = x1;
        x1 = temp;
    }

    if(d1 == 1) // 역이 존재할 경우
        return (x1 > 0 ? x1 : x1+m); // x1이 음수일 경우 x1+m return
    else
        return 0;
}

/*
 * umul_inv() - computes multiplicative inverse a^-1 mod m
 *
 * 입력이 unsigned 64 비트 정수일 때 모듈로 m에서 a의 곱의 역인 a^-1 mod m을 구한다.
 * 만일 역이 존재하지 않으면 0을 리턴한다. 확장유클리드 알고리즘을 변형하여 구현한다.
 * 입출력 모두가 unsigned 64 비트 정수임에 주의한다.
 */
uint64_t umul_inv(uint64_t a, uint64_t m)
{
    uint64_t d0 = a, d1 = m;
    uint64_t x0 = 1, x1 = 0, q, temp;

    while(d1 > 1) {
        q = d0 / d1; // q = d0 div d1

        d0 = d0 - q * d1; // d2 = d0 – q * d1
        temp = d0;
        d0 = d1;
        d1 = temp;

        x0 = x0 - q*x1; // x2 = x0 – q * x1
        temp = x0;
        x0 = x1;
        x1 = temp;
    }

    if(d1 == 1)
        return (x1 >> 63 ? x1+m : x1); // x1이 64비트가 존재할 경우 x1+m, 존재하지 않을 경우 x1 return
    else
        return 0;

}


/*
 * gf16_mul(a, b) - a * b mod x^16+x^5+x^3+x+1
 *
 * 15차식 다항식 a와 b를 곱하고 결과를 16차식 x^16+x^5+x^3+x+1로 나눈 나머지를 계산한다.
 * x^16 = x^5+x^3+x+1 (mod x^16+x^5+x^3+x+1) 특성을 이용한다.
 */

uint16_t xtime(uint16_t x)
{
    // x^15가 존재한다면 x << 1(x^15를 제외한 나머지)와 0x2B(x^5+x^3+x+1를 16진수로 바꿔준 값)를 XOR 한 후 return 
    // 존재하지 않는다면 x << 1을 return 
    return ((x << 1) ^ ((x >> 15) & 1 ? 0x2B : 0));
}

uint16_t gf16_mul(uint16_t a, uint16_t b)
{
    uint16_t r = 0;

    while (b > 0) {
        if (b & 1) r = r ^ a; // b가 1이면 r = r XOR a
        b = b >> 1; // b를 한 비트씩 shift
        a = xtime(a); // xtime 함수 호출
    }
    return r;

}

/*
 * gf16_pow(a,b) - a^b mod x^16+x^5+x^3+x+1
 *
 * 15차식 다항식 a를 b번 지수승한 결과를 16차식 x^16+x^5+x^3+x+1로 나눈 나머지를 계산한다.
 * gf16_mul()과 "Square Multiplication" 알고리즘을 사용하여 구현한다.
 */


uint16_t gf16_pow(uint16_t a, uint16_t b)
{
    uint16_t r = 1;

    while (b > 0) {
        if (b & 1) // b가 1이면 r XOR a를 gf16_mul 함수 이용해 계산
            r = gf16_mul(r, a);
        b = b >> 1; // b를 한 비트씩 shift
        a = gf16_mul(a, a); // gf16_mul 함수 이용해 a값 계산
    }
    return r;

}

/*
 * gf16_inv(a) - a^-1 mod x^16+x^5+x^3+x+1
 *
 * 모둘러 x^16+x^5+x^3+x+1에서 a의 역을 구한다.
 * 역을 구하는 가장 효율적인 방법은 다항식 확장유클리드 알고리즘을 사용하는 것이다.
 * 다만 여기서는 복잡성을 피하기 위해 느리지만 알기 쉬운 지수를 사용하여 구현하였다.
 */
uint16_t gf16_inv(uint16_t a)
{
    return gf16_pow(a, 0xfffe);
}

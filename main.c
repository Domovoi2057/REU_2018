#include <stdint.h>
#include <stdbool.h>
#include <stdio.h>
#include <string.h>
#include "random.h"
#include "fips202.h"
#include "P503_internal.h"
#include "P503_api.h"




extern void multiply32x32_asm (digit_t a, digit_t b, digit_t* c);
extern void multiply256x256_asm(digit_t *a,const digit_t *b,const digit_t *c);


// a : low, b : mid, c : output
extern void test_low_mid(const digit_t *a,const digit_t *b,digit_t *c);
extern void test_mid(const digit_t *a,digit_t *b);
extern void  test_high(const digit_t *a,digit_t *b);
digit_t test_1[16] = { 0x01234567,0x12345678,0x23456789,0x3456789a,0x456789ab,0x56789abc,0x6789abcd,0x789abcde,0x89abcdef,0x9abcdef0,0xabcdef01,0xbcdef012,0xcdef0123,0xdef01234,0xef012345,0xf0123456 };
digit_t test_2[16] = { 0x01234567,0x12345678,0x23456789,0x3456789a,0x456789ab,0x56789abc,0x6789abcd,0x789abcde,0x89abcdef,0x9abcdef0,0xabcdef01,0xbcdef012,0xcdef0123,0xdef01234,0xef012345,0xf0123456 };
digit_t test_3[32];
digit_t test_4[32];
//digit_t test_5[32] ={0x11121314,0x15161718,0x21222324,0x25262728,0x31323334,0x35363738,0x41424344,0x45464748,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF	};
digit_t temp1[16];
digit_t temp2[16];
void multi_testing(const digit_t* a,const digit_t* b,digit_t* c,const unsigned int dummy)
{
	digit_t temp1[16];
digit_t temp2[16];
	// multiply low 
	multiply256x256_asm(temp1,a,b);
	// multiply mid1
	multiply256x256_asm(temp2,a+8,b);
	// add low mid
	test_low_mid(temp1,temp2,c);
	// multiply miid2
	multiply256x256_asm(temp2,b+8,a);
	// add mid 2
	test_mid(temp2,c);
	// multiply high
	multiply256x256_asm(temp2,b+8,a+8);
	// add
	test_high(temp2,c);
}
	
unsigned long cycles(void) { 
       volatile uint32_t *DWT_CYCCNT = (uint32_t *)0xE0001004; #include <stdint.h>
#include <stdbool.h>
#include <stdio.h>
#include <string.h>
#include "random.h"
#include "fips202.h"
#include "P503_internal.h"
#include "P503_api.h"




extern void multiply32x32_asm (digit_t a, digit_t b, digit_t* c);
extern void multiply256x256_asm(digit_t *a,const digit_t *b,const digit_t *c);


// a : low, b : mid, c : output
extern void test_low_mid(const digit_t *a,const digit_t *b,digit_t *c);
extern void test_mid(const digit_t *a,digit_t *b);
extern void  test_high(const digit_t *a,digit_t *b);
digit_t test_1[16] = { 0x01234567,0x12345678,0x23456789,0x3456789a,0x456789ab,0x56789abc,0x6789abcd,0x789abcde,0x89abcdef,0x9abcdef0,0xabcdef01,0xbcdef012,0xcdef0123,0xdef01234,0xef012345,0xf0123456 };
digit_t test_2[16] = { 0x01234567,0x12345678,0x23456789,0x3456789a,0x456789ab,0x56789abc,0x6789abcd,0x789abcde,0x89abcdef,0x9abcdef0,0xabcdef01,0xbcdef012,0xcdef0123,0xdef01234,0xef012345,0xf0123456 };
digit_t test_3[32];
digit_t test_4[32];
//digit_t test_5[32] ={0x11121314,0x15161718,0x21222324,0x25262728,0x31323334,0x35363738,0x41424344,0x45464748,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF	};
digit_t temp1[16];
digit_t temp2[16];
void multi_testing(const digit_t* a,const digit_t* b,digit_t* c,const unsigned int dummy)
{
	digit_t temp1[16];
digit_t temp2[16];
	// multiply low 
	multiply256x256_asm(temp1,a,b);
	// multiply mid1
	multiply256x256_asm(temp2,a+8,b);
	// add low mid
	test_low_mid(temp1,temp2,c);
	// multiply miid2
	multiply256x256_asm(temp2,b+8,a);
	// add mid 2
	test_mid(temp2,c);
	// multiply high
	multiply256x256_asm(temp2,b+8,a+8);
	// add
	test_high(temp2,c);
}
	
unsigned long cycles(void) { 
       volatile uint32_t *DWT_CYCCNT = (uint32_t *)0xE0001004; 
	
        return *DWT_CYCCNT; 
}


int main()
{	
	unsigned char PrivateKeyA[SIDH_SECRETKEYBYTES] = { 0x54, 0x55, 0x55, 0x55,0x55, 0x55, 0x56, 0x78,0x12, 0x34, 0x56, 0x78,0x12, 0x34, 0x56, 0x78,0x12, 0x34, 0x56, 0x78,0x12, 0x34, 0x56, 0x78,0x12, 0x34, 0x56, 0x78,0x12, 0x34, 0x56, 0x78 };//dummy values
	unsigned char PublicKeyA[SIDH_PUBLICKEYBYTES];
	unsigned char PrivateKeyB[SIDH_SECRETKEYBYTES] = { 0x98, 0x76, 0x54, 0x32,0x98, 0x76, 0x54, 0x32,0x98, 0x76, 0x54, 0x32,0x98, 0x76, 0x54, 0x32,0x98, 0x76, 0x54, 0x32,0x98, 0x76, 0x54, 0x32,0x98, 0x76, 0x54, 0x32,0x98, 0x76, 0x54, 0x32 };//dummy values
	unsigned char PublicKeyB[SIDH_PUBLICKEYBYTES];
	unsigned char SharedSecretA[SIDH_BYTES], SharedSecretB[SIDH_BYTES];
		
//	1. Enable the UART module using the RCGCUART register (see page 395).
		
		
//2. Enable the clock to the appropriate GPIO module via the RCGCGPIO register (see page 389).
//To find out which GPIO port to enable, refer to Table 29-5 on page 1932.
//3. Set the GPIO AFSEL bits for the appropriate pins (see page 778). To determine which GPIOs to
//configure, see Table 29-4 on page 1921.
//4. Configure the GPIO current level and/or slew rate as specified for the mode selected (see
//page 780 and page 788).
//5. Configure the PMCn fields in the GPIOPCTL register to assign the UART signals to the appropriate
//pins (see page 795 and Table 29-5 on page 1932).
	
		
		
	EphemeralKeyGeneration_A(PrivateKeyA, PublicKeyA);
	EphemeralKeyGeneration_B(PrivateKeyB, PublicKeyB);
	EphemeralSecretAgreement_A(PrivateKeyA, PublicKeyB, SharedSecretA); 
	EphemeralSecretAgreement_B(PrivateKeyB, PublicKeyA, SharedSecretB);   
	
//	digit_t a1 = 0x01;
//	digit_t b1 = 0xFFFFFFFF;
//	digit_t c1[2] ;
//	digit_t a2 = 0x11111111;
//	digit_t b2 = 0x22222222;
//	digit_t c2[2] ;
//	digit_x_digit(a1,b1,c1);
//	multiply32x32_asm(a2,b2,c2);'
	
	//test_mid(test_1,test_2);
//	multiply256x256_asm(test_3,test_2,test_1);
//	multi_testing(test_1,test_2,test_3,16);
//	mp_mul(test_1,test_2,test_4,16);

	
	
	while (1)
	{
		
	}
	return 1 ;
}
	
        return *DWT_CYCCNT; 
}


int main()
{	
	unsigned char PrivateKeyA[SIDH_SECRETKEYBYTES] = { 0x54, 0x55, 0x55, 0x55,0x55, 0x55, 0x56, 0x78,0x12, 0x34, 0x56, 0x78,0x12, 0x34, 0x56, 0x78,0x12, 0x34, 0x56, 0x78,0x12, 0x34, 0x56, 0x78,0x12, 0x34, 0x56, 0x78,0x12, 0x34, 0x56, 0x78 };//dummy values
	unsigned char PublicKeyA[SIDH_PUBLICKEYBYTES];
	unsigned char PrivateKeyB[SIDH_SECRETKEYBYTES] = { 0x98, 0x76, 0x54, 0x32,0x98, 0x76, 0x54, 0x32,0x98, 0x76, 0x54, 0x32,0x98, 0x76, 0x54, 0x32,0x98, 0x76, 0x54, 0x32,0x98, 0x76, 0x54, 0x32,0x98, 0x76, 0x54, 0x32,0x98, 0x76, 0x54, 0x32 };//dummy values
	unsigned char PublicKeyB[SIDH_PUBLICKEYBYTES];
	unsigned char SharedSecretA[SIDH_BYTES], SharedSecretB[SIDH_BYTES];
		
//	1. Enable the UART module using the RCGCUART register (see page 395).
		
		
//2. Enable the clock to the appropriate GPIO module via the RCGCGPIO register (see page 389).
//To find out which GPIO port to enable, refer to Table 29-5 on page 1932.
//3. Set the GPIO AFSEL bits for the appropriate pins (see page 778). To determine which GPIOs to
//configure, see Table 29-4 on page 1921.
//4. Configure the GPIO current level and/or slew rate as specified for the mode selected (see
//page 780 and page 788).
//5. Configure the PMCn fields in the GPIOPCTL register to assign the UART signals to the appropriate
//pins (see page 795 and Table 29-5 on page 1932).
	
		
		
	EphemeralKeyGeneration_A(PrivateKeyA, PublicKeyA);
	EphemeralKeyGeneration_B(PrivateKeyB, PublicKeyB);
	EphemeralSecretAgreement_A(PrivateKeyA, PublicKeyB, SharedSecretA); 
	EphemeralSecretAgreement_B(PrivateKeyB, PublicKeyA, SharedSecretB);   
	
//	digit_t a1 = 0x01;
//	digit_t b1 = 0xFFFFFFFF;
//	digit_t c1[2] ;
//	digit_t a2 = 0x11111111;
//	digit_t b2 = 0x22222222;
//	digit_t c2[2] ;
//	digit_x_digit(a1,b1,c1);
//	multiply32x32_asm(a2,b2,c2);'
	
	//test_mid(test_1,test_2);
//	multiply256x256_asm(test_3,test_2,test_1);
//	multi_testing(test_1,test_2,test_3,16);
//	mp_mul(test_1,test_2,test_4,16);

	
	
	while (1)
	{
		
	}
	return 1 ;
}

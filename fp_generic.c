/********************************************************************************************
* SIDH: an efficient supersingular isogeny cryptography library
*
* Abstract: portable modular arithmetic for P503
*********************************************************************************************/

#include "P503_internal.h"


// Global constants
extern const uint64_t p503[NWORDS_FIELD];
extern const uint64_t p503p1[NWORDS_FIELD]; 
extern const uint64_t p503x2[NWORDS_FIELD]; 

void fpadd503(const digit_t* a, const digit_t* b, digit_t* c)
{ // Modular addition, c = a+b mod p503.
  // Inputs: a, b in [0, 2*p503-1] 
  // Output: c in [0, 2*p503-1] 
    unsigned int i, carry = 0;
    digit_t mask = 0;
	digit_t temp =0;


	// add
	__asm("adds c[0],a[0],b[0]");
	__asm("adcs c[1],a[1],b[1]");
	__asm("adcs c[2],a[2],b[2]");
	__asm("adcs c[3],a[3],b[3]");
	__asm("adcs c[4],a[4],b[4]");
	__asm("adcs c[5],a[5],b[5]");
	__asm("adcs c[6],a[6],b[6]");
	__asm("adcs c[7],a[7],b[7]");
	__asm("adcs c[8],a[8],b[8]");
	__asm("adcs c[9],a[9],b[9]");
	__asm("adcs c[10],a[10],b[10]");
	__asm("adcs c[11],a[11],b[11]");
	__asm("adcs c[12],a[12],b[12]");
	__asm("adcs c[13],a[13],b[13]");
	__asm("adcs c[14],a[14],b[14]");
	__asm("adcs c[15],a[15],b[15]");
	


		temp = ((digit_t*)p503x2)[0];//subtract p503 from sum.
	__asm("subs c[0],c[0],temp");
		temp = ((digit_t*)p503x2)[1];
	__asm("sbcs c[1],c[1],temp");
		temp = ((digit_t*)p503x2)[2];
	__asm("sbcs c[2],c[2],temp");
		temp = ((digit_t*)p503x2)[3];
	__asm("sbcs c[3],c[3],temp");
		temp = ((digit_t*)p503x2)[4];
	__asm("sbcs c[4],c[4],temp");
		temp = ((digit_t*)p503x2)[5];
	__asm("sbcs c[5],c[5],temp");
		temp = ((digit_t*)p503x2)[6];
	__asm("sbcs c[6],c[6],temp");
		temp = ((digit_t*)p503x2)[7];
	__asm("sbcs c[7],c[7],temp");
		temp = ((digit_t*)p503x2)[8];
	__asm("sbcs c[8],c[8],temp");
		temp = ((digit_t*)p503x2)[9];
	__asm("sbcs c[9],c[9],temp");
		temp = ((digit_t*)p503x2)[10];
	__asm("sbcs c[10],c[10],temp");
		temp = ((digit_t*)p503x2)[11];
	__asm("sbcs c[11],c[11],temp");
		temp = ((digit_t*)p503x2)[12];
	__asm("sbcs c[12],c[12],temp");
		temp = ((digit_t*)p503x2)[13];
	__asm("sbcs c[13],c[13],temp");
		temp = ((digit_t*)p503x2)[14];
	__asm("sbcs c[14],c[14],temp");
		temp = ((digit_t*)p503x2)[15];
	__asm("sbcs c[15],c[15],temp");
	
	
  __asm("sbc mask, #0, #0");//create mask; if carry flag is set, mask  = 0x11


		temp = ((digit_t*)p503x2)[0] & mask;//add (p503AND mask) to difference.
	__asm("adds c[0],c[0],temp");
		temp = ((digit_t*)p503x2)[1] & mask;
	__asm("adcs c[1],c[1],temp");
		temp = ((digit_t*)p503x2)[2] & mask;
	__asm("adcs c[2],c[2],temp");
		temp = ((digit_t*)p503x2)[3] & mask;
	__asm("adcs c[3],c[3],temp");
		temp = ((digit_t*)p503x2)[4] & mask;
	__asm("adcs c[4],c[4],temp");
		temp = ((digit_t*)p503x2)[5] & mask;
	__asm("adcs c[5],c[5],temp");
		temp = ((digit_t*)p503x2)[6] & mask;
	__asm("adcs c[6],c[6],temp");
		temp = ((digit_t*)p503x2)[7] & mask;
	__asm("adcs c[7],c[7],temp");
	temp = ((digit_t*)p503x2)[8] & mask;
	__asm("adcs c[8],c[8],temp");
	temp = ((digit_t*)p503x2)[9] & mask;
	__asm("adcs c[9],c[9],temp");
	temp = ((digit_t*)p503x2)[10] & mask;
	__asm("adcs c[10],c[10],temp");
	temp = ((digit_t*)p503x2)[11] & mask;
	__asm("adcs c[11],c[11],temp");
	temp = ((digit_t*)p503x2)[12] & mask;
	__asm("adcs c[12],c[12],temp");
	temp = ((digit_t*)p503x2)[13] & mask;
	__asm("adcs c[13],c[13],temp");
	temp = ((digit_t*)p503x2)[14] & mask;
	__asm("adcs c[14],c[14],temp");
	temp = ((digit_t*)p503x2)[15] & mask;
	__asm("adcs c[15],c[15],temp");
} 


void fpsub503(const digit_t* a, const digit_t* b, digit_t* c)
{ // Modular subtraction, c = a-b mod p503.
  // Inputs: a, b in [0, 2*p503-1] 
  // Output: c in [0, 2*p503-1] 
    digit_t mask =0, temp =0;

	__asm("subs c[0],a[0],b[0]");
	__asm("sbcs c[1],a[1],b[1]");
	__asm("sbcs c[2],a[2],b[2]");
	__asm("sbcs c[3],a[3],b[3]");
	__asm("sbcs c[4],a[4],b[4]");
	__asm("sbcs c[5],a[5],b[5]");
	__asm("sbcs c[6],a[6],b[6]");
	__asm("sbcs c[7],a[7],b[7]");
	__asm("sbcs c[8],a[8],b[8]");
	__asm("sbcs c[9],a[9],b[9]");
	__asm("sbcs c[10],a[10],b[10]");
	__asm("sbcs c[11],a[11],b[11]");
	__asm("sbcs c[12],a[12],b[12]");
	__asm("sbcs c[13],a[13],b[13]");
	__asm("sbcs c[14],a[14],b[14]");
	__asm("sbcs c[15],a[15],b[15]");//subtract 16 times
	
	
	__asm("sbc mask, #0, #0");//create mask; if carry flag is set, mask  = 0x11
			
	
	temp = ((digit_t*)p503x2)[0] & mask;//add (2*p503 AND mask) to difference.
	__asm("adds c[0],c[0],temp");
		temp = ((digit_t*)p503x2)[1] & mask;
	__asm("adcs c[1],c[1],temp");
		temp = ((digit_t*)p503x2)[2] & mask;
	__asm("adcs c[2],c[2],temp");
		temp = ((digit_t*)p503x2)[3] & mask;
	__asm("adcs c[3],c[3],temp");
		temp = ((digit_t*)p503x2)[4] & mask;
	__asm("adcs c[4],c[4],temp");
		temp = ((digit_t*)p503x2)[5] & mask;
	__asm("adcs c[5],c[5],temp");
		temp = ((digit_t*)p503x2)[6] & mask;
	__asm("adcs c[6],c[6],temp");
		temp = ((digit_t*)p503x2)[7] & mask;
	__asm("adcs c[7],c[7],temp");
	temp = ((digit_t*)p503x2)[8] & mask;
	__asm("adcs c[8],c[8],temp");
	temp = ((digit_t*)p503x2)[9] & mask;
	__asm("adcs c[9],c[9],temp");
	temp = ((digit_t*)p503x2)[10] & mask;
	__asm("adcs c[10],c[10],temp");
	temp = ((digit_t*)p503x2)[11] & mask;
	__asm("adcs c[11],c[11],temp");
	temp = ((digit_t*)p503x2)[12] & mask;
	__asm("adcs c[12],c[12],temp");
	temp = ((digit_t*)p503x2)[13] & mask;
	__asm("adcs c[13],c[13],temp");
	temp = ((digit_t*)p503x2)[14] & mask;
	__asm("adcs c[14],c[14],temp");
	temp = ((digit_t*)p503x2)[15] & mask;
	__asm("adcs c[15],c[15],temp");
} 


void fpsub503_original(const digit_t* a, const digit_t* b, digit_t* c)
{ // Modular subtraction, c = a-b mod p503.
  // Inputs: a, b in [0, 2*p503-1] 
  // Output: c in [0, 2*p503-1] 
    unsigned int i, borrow = 0;
    digit_t mask;

    for (i = 0; i < NWORDS_FIELD; i++) {
        SUBC(borrow, a[i], b[i], borrow, c[i]); 
    }
    mask = 0 - (digit_t)borrow;

    borrow = 0;
    for (i = 0; i < NWORDS_FIELD; i++) {
        ADDC(borrow, c[i], ((digit_t*)p503x2)[i] & mask, borrow, c[i]); 
    }
}


void fpneg503(digit_t* a)
{ // Modular negation, a = -a mod p503.
  // Input/output: a in [0, 2*p503-1] 
    unsigned int i, borrow = 0;

    for (i = 0; i < NWORDS_FIELD; i++) {
        SUBC(borrow, ((digit_t*)p503x2)[i], a[i], borrow, a[i]); 
    }
}


void fpdiv2_503(const digit_t* a, digit_t* c)
{ // Modular division by two, c = a/2 mod p503.
  // Input : a in [0, 2*p503-1] 
  // Output: c in [0, 2*p503-1] 
    unsigned int i, carry = 0;
    digit_t mask;
        
    mask = 0 - (digit_t)(a[0] & 1);    // If a is odd compute a+p503
    for (i = 0; i < NWORDS_FIELD; i++) {
        ADDC(carry, a[i], ((digit_t*)p503)[i] & mask, carry, c[i]); 
    }

    mp_shiftr1(c, NWORDS_FIELD);
} 


void fpcorrection503(digit_t* a)
{ // Modular correction to reduce field element a in [0, 2*p503-1] to [0, p503-1].
    unsigned int i, borrow = 0;
    digit_t mask;

    for (i = 0; i < NWORDS_FIELD; i++) {
        SUBC(borrow, a[i], ((digit_t*)p503)[i], borrow, a[i]); 
    }
    mask = 0 - (digit_t)borrow;

    borrow = 0;
    for (i = 0; i < NWORDS_FIELD; i++) {
        ADDC(borrow, a[i], ((digit_t*)p503)[i] & mask, borrow, a[i]); 
    }
}

/*
void digit_x_digit(const digit_t a, const digit_t b, digit_t* c)
{ // Digit multiplication, digit * digit -> 2-digit result    
    register digit_t al, ah, bl, bh, temp;
    digit_t albl, albh, ahbl, ahbh, res1, res2, res3, carry;
    digit_t mask_low = (digit_t)(-1) >> (sizeof(digit_t)*4), mask_high = (digit_t)(-1) << (sizeof(digit_t)*4);

    al = a & mask_low;                        // Low part
    ah = a >> (sizeof(digit_t) * 4);          // High part
    bl = b & mask_low;
    bh = b >> (sizeof(digit_t) * 4);

    albl = al*bl;
    albh = al*bh;
    ahbl = ah*bl;
    ahbh = ah*bh;
    c[0] = albl & mask_low;                   // C00

    res1 = albl >> (sizeof(digit_t) * 4);
    res2 = ahbl & mask_low;
    res3 = albh & mask_low;  
    temp = res1 + res2 + res3;
    carry = temp >> (sizeof(digit_t) * 4);
    c[0] ^= temp << (sizeof(digit_t) * 4);    // C01   

    res1 = ahbl >> (sizeof(digit_t) * 4);
    res2 = albh >> (sizeof(digit_t) * 4);
    res3 = ahbh & mask_low;
    temp = res1 + res2 + res3 + carry;
    c[1] = temp & mask_low;                   // C10 
    carry = temp & mask_high; 
    c[1] ^= (ahbh & mask_high) + carry;       // C11
}*/


//extern void multiply32x32_asm (digit_t a, digit_t b, digit_t* c);
//extern void multiply32x32_asm_m4 (digit_t a, digit_t b, digit_t* c);
#define digit_x_digit multiply32x32_asm_m4



void mp_mul(const digit_t* a, const digit_t* b, digit_t* c, const unsigned int nwords)//not used in program.
{ // Multiprecision comba multiply, c = a*b, where lng(a) = lng(b) = nwords.   
    unsigned int i, j;
    digit_t t = 0, u = 0, v = 0, UV[2];
    unsigned int carry = 0;
    
    for (i = 0; i < nwords; i++) {
        for (j = 0; j <= i; j++) {
            MUL(a[j], b[i-j], UV+1, UV[0]); 
            ADDC(0, UV[0], v, carry, v); 
            ADDC(carry, UV[1], u, carry, u); 
            t += carry;
        }
        c[i] = v;
        v = u; 
        u = t;
        t = 0;
    }

    for (i = nwords; i < 2*nwords-1; i++) {
        for (j = i-nwords+1; j < nwords; j++) {
            MUL(a[j], b[i-j], UV+1, UV[0]); 
            ADDC(0, UV[0], v, carry, v); 
            ADDC(carry, UV[1], u, carry, u); 
            t += carry;
        }
        c[i] = v;
        v = u; 
        u = t;
        t = 0;
    }
    c[2*nwords-1] = v; 
}

//extern void test_loop(digit_t *m);

void rdc_mont(const digit_t* ma, digit_t* mc)
{ // Efficient Montgomery reduction using comba and exploiting the special form of the prime p503.
  // mc = ma*R^-1 mod p503x2, where R = 2^512.
  // If ma < 2^512*p503, the output mc is in the range [0, 2*p503-1].
  // ma is assumed to be in Montgomery representation.
    unsigned int i, j, carry, count = p503_ZERO_WORDS;
    digit_t UV[2], t = 0, u = 0, v = 0;
		digit_t temp1;
    //for (i = 0; i < NWORDS_FIELD; i++) {
    //    mc[i] = 0;
    //}
	//test_loop(mc);	
	
    for (i = 0; i < NWORDS_FIELD; i++) {
        for (j = 0; j < i; j++) {
            if (j < (i-p503_ZERO_WORDS+1)) { 
							temp1 = ((digit_t*)p503p1)[i-j];
               // MUL(mc[j], ((digit_t*)p503p1)[i-j], UV+1, UV[0]);
								__asm("umull UV[0],UV[1],mc[j],temp1");
                //ADDC(0, UV[0], v, carry, v); 
                //ADDC(carry, UV[1], u, carry, u); 
                //t += carry; 
							__asm("adds v,v,UV[0]");
							__asm("adcs u,u,UV[1]");
							__asm("adcs t,t,#0x00");
            }
        }
        //ADDC(0, v, ma[i], carry, v); 
        //ADDC(carry, u, 0, carry, u); 				
       //t += carry; 
				__asm("adds v,v,ma[i]");
				__asm("adcs u,u,#0x00");
				__asm("adcs t,t,#0x00");
        mc[i] = v;
        v = u;
        u = t;
        t = 0;
    }    

    for (i = NWORDS_FIELD; i < 2*NWORDS_FIELD-1; i++) {
        if (count > 0) {
            count -= 1;
        }
        for (j = i-NWORDS_FIELD+1; j < NWORDS_FIELD; j++) {
            if (j < (NWORDS_FIELD-count)) { 
							temp1 = ((digit_t*)p503p1)[i-j];
              //  MUL(mc[j], ((digit_t*)p503p1)[i-j], UV+1, UV[0]);
							__asm("umull UV[0],UV[1],mc[j],temp1");
                //ADDC(0, UV[0], v, carry, v); 
                //ADDC(carry, UV[1], u, carry, u); 
                //t += carry;
							__asm("adds v,v,UV[0]");
							__asm("adcs u,u,UV[1]");
							__asm("adcs t,t,#0x00");
            }
        }
        //ADDC(0, v, ma[i], carry, v); 
        //ADDC(carry, u, 0, carry, u); 
        //t += carry; 
				__asm("adds v,v,ma[i]");
				__asm("adcs u,u,#0x00");
				__asm("adcs t,t,#0x00");
        mc[i-NWORDS_FIELD] = v;
        v = u;
        u = t;
        t = 0;
    }
    //ADDC(0, v, ma[2*NWORDS_FIELD-1], carry, v); 
		__asm("adds v,v,ma[31]");		// ma[2*NWORDS_FIELD-1]
    //mc[NWORDS_FIELD-1] = v;
		mc[15] = v;
}
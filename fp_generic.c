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

void rdc_mont1(const digit_t* ma, digit_t* mc)
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
	
    for (i = 0; i < NWORDS_FIELD; i++) {		//round 1, description:
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

    for (i = NWORDS_FIELD; i < 2*NWORDS_FIELD-1; i++) {			//round 2, description:
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

void rdc_mont(const digit_t* ma, digit_t* mc)
{ // Efficient Montgomery reduction using comba and exploiting the special form of the prime p503.
//  // mc = ma*R^-1 mod p503x2, where R = 2^512.
//  // If ma < 2^512*p503, the output mc is in the range [0, 2*p503-1].
//  // ma is assumed to be in Montgomery representation.
//    unsigned int i, j, carry, count = p503_ZERO_WORDS;
    digit_t UV[2], t = 0, u = 0, v = 0;
		digit_t temp1;
//    //for (i = 0; i < NWORDS_FIELD; i++) {
//    //    mc[i] = 0;
//    //}
//	//test_loop(mc);	
		//round 1
		__asm("mov mc[0], ma[0]");//	round 1: ii=0, jj=0; there shouldn't be any initial carries, so the propagation is unnecessary and entire iteration can be written as a move op.
		__asm("mov mc[1], ma[1]");//	round 1: ii=1, jj=0; there shouldn't be any initial carries, so the propagation is unnecessary and entire iteration can be written as a move op.
		__asm("mov mc[2], ma[2]");//	round 1: ii=2, jj=1; there shouldn't be any initial carries, so the propagation is unnecessary and entire iteration can be written as a move op.
		__asm("mov mc[3], ma[3]");//	round 1: ii=3, jj=2; there shouldn't be any initial carries, so the propagation is unnecessary and entire iteration can be written as a move op.
		__asm("mov mc[4], ma[4]");//	round 1: ii=4, jj=3; there shouldn't be any initial carries, so the propagation is unnecessary and entire iteration can be written as a move op.
		__asm("mov mc[5], ma[5]");//	round 1: ii=5, jj=4; there shouldn't be any initial carries, so the propagation is unnecessary and entire iteration can be written as a move op.
		__asm("mov mc[6], ma[6]");//	round 1: ii=6, jj=5; there shouldn't be any initial carries, so the propagation is unnecessary and entire iteration can be written as a move op.
		__asm("mov mc[7], ma[7]");//	round 1: ii=7, jj=6; there shouldn't be any initial carries, so the propagation is unnecessary and entire iteration can be written as a move op.
		__asm("mov mc[8], ma[8]");//	round 1: ii=8, jj=7; there shouldn't be any initial carries, so the propagation is unnecessary and entire iteration can be written as a move op.

		temp1 = ((digit_t*)p503p1)[9];
		__asm("umull UV[0], UV[1], temp1, mc[0]");//	round 1: ii=9, jj=0; multiply 10th dgt of P503p1 w/ 1st dgt of C.
		__asm("adds mc[9], UV[0], ma[9]");//	round 1: ii=9, jj=9; add LSB w/ 10th dgt of A. Store in 10th dgt of C. Set carry.
		__asm("adcs v, UV[1], #0x00");//	add MSB w/ 0, store in v. Set carry.
		__asm("adcs u, #0x00, #0x00");//	Add 0 to 0, to set carry flag. Store in u.
		
		temp1 = ((digit_t*)p503p1)[10];		
		__asm("umull UV[0], UV[1], temp1, mc[0]");// round 1: ii=10, jj=0; multiply 11th dgt of P503p1 w/ 1st dgt of C.
		__asm("adds v, v, UV[0]");// add LSB to v, store in v.
		__asm("adcs u, u, UV[1]");// add MSB to u, store in u.
		__asm("adcs t, #0x00, #0x00");//	place carry in t.

		temp1 = ((digit_t*)p503p1)[9];
		__asm("umull UV[0], UV[1], temp1, mc[1]");//	round 1: ii=10, jj=1; multiply 10th dgt of P503p1 w/ 2nd dgt of C.
		__asm("adds v, v, UV[0]");//	add LSB to v, store in v.
		__asm("adcs u, u, UV[1]");//	add MSB to u, store in u.
		__asm("adcs t, t, #0x00");//	add carry to t.
		
		__asm("adds mc[10], v, ma[10]");//	round 1: ii=10, jj=10; add 11th dgt of A to v, store in 11th dgt of C.
		__asm("adcs v, u, #0x00");//	add carry to u, store sum in v. Propagate carry.
		__asm("adcs u, t, #0x00");//	add carry to t, store sum in u.
		__asm("mov t, #0x00");//	set t to zero.
		
		
		temp1 = ((digit_t*)p503p1)[11];		
		__asm("umull UV[0], UV[1], temp1, mc[0]");//	round 1: ii=11, jj=0; multiply 12th dgt of P503p1 w/ 1st dgt of C.
		__asm("adds v, v, UV[0]");//	add LSB to v, store in v.
		__asm("adcs u, u, UV[1]");//	add MSB to u, store in u.
		__asm("adcs t, t, #0x00");//	add carry to t.

		temp1 = ((digit_t*)p503p1)[10];
		__asm("umull UV[0], UV[1], temp1, mc[1]");//	round 1: ii=11, jj=1; multiply 11th digit of P503p1 w/ 2nd digit of C.
		__asm("adds v, v, UV[0]");//	add LSB to v, store in v.
		__asm("adcs u, u, UV[1]");// add MSB to u, store in u.
		__asm("adcs t, t, #0x00");//	add carry to t.

		temp1 = ((digit_t*)p503p1)[9];
		__asm("umull UV[0], UV[1], temp1, mc[2]");//	round 1: ii=11, jj=2; multiply 10th dgt of p503p1 w/ 3rd dgt of C.
		__asm("adds v, v, UV[0]");//	add LSB to v, store in v.
		__asm("adcs u, u, UV[1]");//	add MSB to u, store in u.
		__asm("adcs t, t, #0x00");// add carry to t.
		
		__asm("adds mc[11], v, ma[11]");//	round 1: ii=11, jj=11; add v to 12th dgt of A, store in 12th dgt of C.
		__asm("adcs v, u, #0x00");//	add carry to u, store sum in v. Propagate carry.
		__asm("adcs u, t, #0x00");// add carry to t, store sum in u. Propagate carry.
		__asm("mov t, #0x00");//	set t to zero.
		

		temp1 = ((digit_t*)p503p1)[12];
		__asm("umull UV[0], UV[1], temp1, mc[0]");//	round 1: ii=12, jj=0; multiply 13th dgt of p503p1 w/ 1st dgt of C.
		__asm("adds v, v, UV[0]");//	add LSB to v, store in v.
		__asm("adcs u, u, UV[1]");//	add MSB to u, store in u.
		__asm("adcs t, t, #0x00");//	add carry to t, store in t.

		temp1 = ((digit_t*)p503p1)[11];
		__asm("umull UV[0], UV[1], temp1, mc[1]");//	round 1: ii=12, jj=1; multiply 12th dgt of p503p1 w/ 2nd dgt of C.
		__asm("adds v, v, UV[0]");//	add LSB to v, store in v.
		__asm("adcs u, u, UV[1]");//	add MSB to u, store in u.
		__asm("adcs t, t, #0x00");//	add carry to t, store in t.

		temp1 = ((digit_t*)p503p1)[10];
		__asm("umull UV[0], UV[1], temp1, mc[2]");//	round 1: ii=12, jj=2; multiply 11th dgt of p503p1 w/ 3rd dgt of C.
		__asm("adds v, v, UV[0]");//	add LSB to v, store in v.
		__asm("adcs u, u, UV[1]");//	add MSB to u, store in u.
		__asm("adcs t, t, #0x00");//	add carry to t, store in t.

		temp1 = ((digit_t*)p503p1)[9];
		__asm("umull UV[0], UV[1], temp1, mc[3]");//	round 1: ii=12, jj=3; multiply 10th dgt of p503p1 w/ 4th dgt of C.
		__asm("adds v, v, UV[0]");//	add LSB to v, store in v.
		__asm("adcs u, u, UV[1]");//	add MSB to u, store in u.
		__asm("adcs t, t, #0x00");//	add carry to t, store in t.
		
		__asm("adds mc[12], v, ma[12]");//	round 1: ii=12, jj=12; add 13th dgt of A to v, store in 13th dgt of C.
		__asm("adcs v, u, #0x00");//	add carry to u, store sum in v.
		__asm("adcs u, t, #0x00");// add carry to t, store in u.
		__asm("mov t, #0x00");//	set t to zero.
		
		
		temp1 = ((digit_t*)p503p1)[13];		
		__asm("umull UV[0], UV[1], temp1, mc[0]");//	round 1: ii=13, jj=0;	multiply 14th dgt of p503p1 w/ 1st dgt of C.
		__asm("adds v, v, UV[0]");//	add LSB to v, store in v.
		__asm("adcs u, u, UV[1]");//	add MSB to u, store in u.
		__asm("adcs t, t, #0x00");//	add carry to t, store in t.
	
		temp1 = ((digit_t*)p503p1)[12];	
		__asm("umull UV[0], UV[1], temp1, mc[1]");//	round 1: ii=13, jj=1; multiply 13th dgt of p503p1 w/ 2nd dgt of C.
		__asm("adds v, v, UV[0]");//	add LSB to v, store in v.
		__asm("adcs u, u, UV[1]");//	add MSB to u, store in u.
		__asm("adcs t, t, #0x00");//	add carry to t, store in t.
		
		temp1 = ((digit_t*)p503p1)[11];
		__asm("umull UV[0], UV[1], temp1, mc[2]");//	round 1: ii=13, jj=2; multiply 12th dgt of p503p1 w/ 3rd dgt of C.
		__asm("adds v, v, UV[0]");//	add LSB to v, store in v.
		__asm("adcs u, u, UV[1]");//	add MSB to u, store in u.
		__asm("adcs t, t, #0x00");//	add carry to t, store in t.
		
		temp1 = ((digit_t*)p503p1)[10];		
		__asm("umull UV[0], UV[1], temp1, mc[3]");//	round 1: ii=13, jj=3; multiply 11th dgt of p503p1 w/ 4th dgt of C.
		__asm("adds v, v, UV[0]");//	add LSB to v, store in v.
		__asm("adcs u, u, UV[1]");//	add MSB to u, store in u.
		__asm("adcs t, t, #0x00");//	add carry to t, store in t.

		temp1 = ((digit_t*)p503p1)[9];
		__asm("umull UV[0], UV[1], temp1, mc[4]");//	round 1: ii=13, jj=4; multiply 10th dgt of p503p1 w/ 5th dgt of C.
		__asm("adds v, v, UV[0]");//	add LSB to v, store in v.
		__asm("adcs u, u, UV[1]");//	add MSB to u, store in u.
		__asm("adcs t, t, #0x00");//	add carry to t, store in t.
		
		__asm("adds mc[13], v, ma[13]");//	round 1: ii=13, jj=13; add v to 14th dgt of A, store in 14th dgt of C.
		__asm("adcs v, u, #0x00");//	add carry to u, store sum in v.
		__asm("adcs u, t, #0x00");//	add carry to t, store sum in u.
		__asm("mov t, #0x00");// set t to zero.
		
		
		temp1 = ((digit_t*)p503p1)[14];		
		__asm("umull UV[0], UV[1], temp1, mc[0]");//	round 1: ii=14, jj=0; multiply 15th dgt of p503p1 w/ 1st dgt of C.
		__asm("adds v, v, UV[0]");//	add LSB to v, store in v.
		__asm("adcs u, u, UV[1]");//	add MSB to u, store in u.
		__asm("adcs t, t, #0x00");// add carry to t, store in t.

		temp1 = ((digit_t*)p503p1)[13];
		__asm("umull UV[0], UV[1], temp1, mc[1]");//	round 1: ii=14, jj=1; multiply 14th dgt of p503p1 w/ 2nd dgt of C.
		__asm("adds v, v, UV[0]");//	add LSB to v, store in v.
		__asm("adcs u, u, UV[1]");//	add MSB to u, store in u.
		__asm("adcs t, t, #0x00");//	add carry to t, store in t.

		temp1 = ((digit_t*)p503p1)[12];
		__asm("umull UV[0], UV[1], temp1, mc[2]");//	round 1: ii=14, jj=2; multiply 13th dgt of p503p1 w/ 3rd dgt of C.
		__asm("adds v, v, UV[0]");//	add LSB to v, store in v.
		__asm("adcs u, u, UV[1]");// add MSB to u, store in u.
		__asm("adcs t, t, #0x00");//	 add carry to t, store in t.

		temp1 = ((digit_t*)p503p1)[11];
		__asm("umull UV[0], UV[1], temp1, mc[3]");//	round 1: ii=14, jj=3; multiply 12th dgt of p503p1 w/ 4th dgt of C.
		__asm("adds v, v, UV[0]");// add LSB to v, store in v.
		__asm("adcs u, u, UV[1]");//	add MSB to u, store in u.
		__asm("adcs t, t, #0x00");//	add carry to t, store in t.

		temp1 = ((digit_t*)p503p1)[10];
		__asm("umull UV[0], UV[1], temp1, mc[4]");//	round 1: ii=14, jj=4; multiply 11th dgt of p503p1 w/ 5th dgt of C.
		__asm("adds v, v, UV[0]");// add LSB to v, store in v.
		__asm("adcs u, u, UV[1]");// add MSB to u, store in u.
		__asm("adcs t, t, #0x00");//	add carry to t, store in t.

		temp1 = ((digit_t*)p503p1)[9];
		__asm("umull UV[0], UV[1], temp1, mc[5]");// round 1: ii=14, jj=5; muiltiply 10th dgt of p503p1 w/ 6th dgt of C.
		__asm("adds v, v, UV[0]");//	add LSB to v, store in v.
		__asm("adcs u, u, UV[1]");//	add MSB to u, store in u.
		__asm("adcs t, t, #0x00");//	add carry to t, store in t.
		
		__asm("adds mc[14], v, ma[14]");//	round 1: ii=14, jj=14; add 15th dgt of A to v, store in 15th dgt of C.
		__asm("adcs v, u, #0x00");//	add carry to u, store in v.
		__asm("adcs u, t, #0x00");//	add carry to t, store in u.
		__asm("mov t, #0x00");//	set t to zero.
		
		
		temp1 = ((digit_t*)p503p1)[15];		
		__asm("umull UV[0], UV[1], temp1, mc[0]");//	round 1: ii=15, jj=0; multiply 16th dgt of p503p1 w/ 1st dgt of C.
		__asm("adds v, v, UV[0]");//	add LSB to v, store in v.
		__asm("adcs u, u, UV[1]");//	add MSB to u, store in u.
		__asm("adcs t, t, #0x00");//	add carry to t, store in t.

		temp1 = ((digit_t*)p503p1)[14];
		__asm("umull UV[0], UV[1], temp1, mc[1]");//	round 1: ii=15, jj=1; multiply 15th dgt of p503p1 w/ 2nd dgt of C.
		__asm("adds v, v, UV[0]");//	add LSB to v, store in v.
		__asm("adcs u, u, UV[1]");//	add MSB to u, store in u.
		__asm("adcs t, t, #0x00");//	add carry to t, store in t.

		temp1 = ((digit_t*)p503p1)[13];
		__asm("umull UV[0], UV[1], temp1, mc[2]");//	round 1: ii=15, jj=2; multiply 14th dgt of p503p1 w/ 3rd dgt of C.
		__asm("adds v, v, UV[0]");//	add LSB to v, store in v.
		__asm("adcs u, u, UV[1]");//	add MSB to u, store in u.
		__asm("adcs t, t, #0x00");//	add carry to t, store in t.
		
		temp1 = ((digit_t*)p503p1)[12];		
		__asm("umull UV[0], UV[1], temp1, mc[3]");//	round 1: ii=15, jj=3; multiply 13th dgt of p503p1 w/ 4th dgt of C.
		__asm("adds v, v, UV[0]");//	add LSB to v, store in v.
		__asm("adcs u, u, UV[1]");// add MSB to u, store in u.
		__asm("adcs t, t, #0x00");//	add carry to t, store in t.

		temp1 = ((digit_t*)p503p1)[11];
		__asm("umull UV[0], UV[1], temp1, mc[4]");//	round 1: ii=15, jj=4; multiply 12th dgt of p503p1 w/ 5th dgt of C.
		__asm("adds v, v, UV[0]");//	add LSB to v, store in v.
		__asm("adcs u, u, UV[1]");//	add MSB to u, store in u.
		__asm("adcs t, t, #0x00");//	add carry to t, store in t.

		temp1 = ((digit_t*)p503p1)[10];
		__asm("umull UV[0], UV[1], temp1, mc[5]");//	round 1: ii=15, jj=5; multiply 11th dgt of p503p1 w/ 6th dgt of C.
		__asm("adds v, v, UV[0]");//	add LSB to v, store in v.
		__asm("adcs u, u, UV[1]");//	add MSB to u, store in u.
		__asm("adcs t, t, #0x00");//	add carry to t, store in t.
		
		temp1 = ((digit_t*)p503p1)[9];		
		__asm("umull UV[0], UV[1], temp1, mc[6]");//	round 1: ii=15, jj=6; multiply 10th dgt of p503p1 w/ 7th dgt of C.
		__asm("adds v, v, UV[0]");//	add LSB to v, store in v.
		__asm("adcs u, u, UV[1]");//	add MSB to u, store in u.
		__asm("adcs t, t, #0x00");//	add carry to t, store in t.
		
		__asm("adds mc[15], v, ma[15]");//	round 1: ii=15, jj=15; add 16th dgt of A to v, store in 16th dgt of C.
		__asm("adcs v, u, #0x00");//	add carry to u, store in v.
		__asm("adcs u, t, #0x00");//	add carry to t, store in u.
		__asm("mov t, #0x00");//	set t to zero.
		
		
		temp1 = ((digit_t*)p503p1)[15];		
		__asm("umull UV[0], UV[1], temp1, mc[1]");//	round 2: ii=16, jj=1; multiply 16th dgt of p503p1 w/ 2nd dgt of C.
		__asm("adds v, v, UV[0]");//	add LSB to v, store in v.
		__asm("adcs u, u, UV[1]");//	add MSB to u, store in u.
		__asm("adcs t, t, #0x00");//	add carry to t, store in t.

		temp1 = ((digit_t*)p503p1)[14];
		__asm("umull UV[0], UV[1], temp1, mc[2]");//	round 2: ii=16, jj=2; multiply 15th dgt of p503p1 w/ 3rd dgt of C.
		__asm("adds v, v, UV[0]");//	add LSB to v, store in v.
		__asm("adcs u, u, UV[1]");//	add MSB to u, store in u.
		__asm("adcs t, t, #0x00");//	add carry to t, store in t.

		temp1 = ((digit_t*)p503p1)[13];
		__asm("umull UV[0], UV[1], temp1, mc[3]");//	round 2: ii=16, jj=3; multiply 14th dgt of p503p1 w/ 4th dgt of C.
		__asm("adds v, v, UV[0]");//	add LSB to v, store in v.
		__asm("adcs u, u, UV[1]");//	add MSB to u, store in u.
		__asm("adcs t, t, #0x00");//	add carry to t, store in t.

		temp1 = ((digit_t*)p503p1)[12];
		__asm("umull UV[0], UV[1], temp1, mc[4]");//	round 2: ii=16, jj=4; multiply 13th dgt of p503p1 w/ 5th dgt of C.
		__asm("adds v, v, UV[0]");//	add LSB to v, store in v.
		__asm("adcs u, u, UV[1]");//	add MSB to u, store in u.
		__asm("adcs t, t, #0x00");//	add carry to t, store in t.

		temp1 = ((digit_t*)p503p1)[11];
		__asm("umull UV[0], UV[1], temp1, mc[5]");//	round 2: ii=16, jj=5; multiply 12th dgt of p503p1 w/ 6th dgt of C.
		__asm("adds v, v, UV[0]");//	add LSB to v, store in v.
		__asm("adcs u, u, UV[1]");//	add MSB to u, store in u.
		__asm("adcs t, t, #0x00");//	add carry to t, store in t.

		temp1 = ((digit_t*)p503p1)[10];
		__asm("umull UV[0], UV[1], temp1, mc[6]");//	round 2: ii=16, jj=6; multiply 11th dgt of p503p1 w/ 7th dgt of C.
		__asm("adds v, v, UV[0]");//	add LSB to v, store in v.
		__asm("adcs u, u, UV[1]");//	add MSB to u, store in u.
		__asm("adcs t, t, #0x00");//	add carry to t, store in t.

		temp1 = ((digit_t*)p503p1)[9];
		__asm("umull UV[0], UV[1], temp1, mc[7]");//	round 2: ii=16, jj=7; multiply 10th dgt of p503p1 w/ 8th dgt of C.
		__asm("adds v, v, UV[0]");//	add LSB to v, store in v.
		__asm("adcs u, u, UV[1]");//	add MSB to u, store in u.
		__asm("adcs t, t, #0x00");//	add carry to t, store in t.
	
		temp1 = ((digit_t*)p503p1)[8];	
		__asm("umull UV[0], UV[1], temp1, mc[8]");//	round 2: ii=16, jj=8; multiply 9th dgt of p503p1 w/ 9th dgt of C.
		__asm("adds v, v, UV[0]");//	add LSB to v, store in v.
		__asm("adcs u, u, UV[1]");//	add MSB to u, store in u.
		__asm("adcs t, t, #0x00");//	add carry to t, store in t.

		temp1 = ((digit_t*)p503p1)[7];
		__asm("umull UV[0], UV[1], temp1, mc[9]");//	round 2: ii=16, jj=9; multiply 8th dgt of p503p1 w/ 10th dgt of C.
		__asm("adds v, v, UV[0]");//	add LSB to v, store in v.
		__asm("adcs u, u, UV[1]");//	add MSB to u, store in u.
		__asm("adcs t, t, #0x00");//	add carry to t, store in t.

		__asm("adds mc[0], v, ma[16]");//	round 2: ii=16, jj=16;;	add v to 17th dgt of A, store in 1st dgt of C.
		__asm("adcs v, u, #0x00");//	add carry to u, store in v.
		__asm("adcs u, t, #0x00");//	add carry to t, store in u.
		__asm("mov t, #0x00");//	set t to zero.
		
		
		temp1 = ((digit_t*)p503p1)[15];		
		__asm("umull UV[0], UV[1], temp1, mc[2]");//	round 2: ii=17, jj=2;	multiply 16th dgt of p503p1 w/ 3rd dgt of C.
		__asm("adds v, v, UV[0]");//	add LSB to v, store in v.
		__asm("adcs u, u, UV[1]");//	add MSB to u, store in u.
		__asm("adcs t, t, #0x00");//	add carry to t, store in t.
		
		temp1 = ((digit_t*)p503p1)[14];		
		__asm("umull UV[0], UV[1], temp1, mc[3]");//	round 2: ii=17, jj=3;	multiply 15th dgt of p503p1 w/ 4th dgt of C.
		__asm("adds v, v, UV[0]");//	add LSB to v, store in v.
		__asm("adcs u, u, UV[1]");//	add MSB to u, store in u.
		__asm("adcs t, t, #0x00");//	add carry to t, store in t.

		temp1 = ((digit_t*)p503p1)[13];
		__asm("umull UV[0], UV[1], temp1, mc[4]");//	round 2: ii=17, jj=4;	multiply 14th dgt of p503p1 w/ 5th dgt of C.
		__asm("adds v, v, UV[0]");//	add LSB to v, store in v.
		__asm("adcs u, u, UV[1]");//	add MSB to u, store in u.
		__asm("adcs t, t, #0x00");//	add carry to t, store in t.

		temp1 = ((digit_t*)p503p1)[12];
		__asm("umull UV[0], UV[1], temp1, mc[5]");//	round 2: ii=17, jj=5;	multiply 13th dgt of p503p1 w/ 6th dgt of C.
		__asm("adds v, v, UV[0]");//	add LSB to v, store in v.
		__asm("adcs u, u, UV[1]");//	add MSB to u, store in u.
		__asm("adcs t, t, #0x00");//	add carry to t, store in t.		

		temp1 = ((digit_t*)p503p1)[11];
		__asm("umull UV[0], UV[1], temp1, mc[6]");//	round 2: ii=17, jj=6;	multiply 12th dgt of p503p1 w/ 7th dgt of C.
		__asm("adds v, v, UV[0]");//	add LSB to v, store in v.
		__asm("adcs u, u, UV[1]");//	add MSB to u, store in u.
		__asm("adcs t, t, #0x00");//	add carry to t, store in t.		
		
		temp1 = ((digit_t*)p503p1)[10];		
		__asm("umull UV[0], UV[1], temp1, mc[7]");//	round 2: ii=17, jj=7;	multiply 11th dgt of p503p1 w/ 8th dgt of C.
		__asm("adds v, v, UV[0]");//	add LSB to v, store in v.
		__asm("adcs u, u, UV[1]");//	add MSB to u, store in u.
		__asm("adcs t, t, #0x00");//	add carry to t, store in t.		
		
		temp1 = ((digit_t*)p503p1)[9];		
		__asm("umull UV[0], UV[1], temp1, mc[8]");//	round 2: ii=17, jj=8;	multiply 10th dgt of p503p1 w/ 9th dgt of C.
		__asm("adds v, v, UV[0]");//	add LSB to v, store in v.
		__asm("adcs u, u, UV[1]");//	add MSB to u, store in u.
		__asm("adcs t, t, #0x00");//	add carry to t, store in t.		

		temp1 = ((digit_t*)p503p1)[8];
		__asm("umull UV[0], UV[1], temp1, mc[9]");//	round 2: ii=17, jj=9;	multiply 9th dgt of p503p1 w/ 10th dgt of C.
		__asm("adds v, v, UV[0]");//	add LSB to v, store in v.
		__asm("adcs u, u, UV[1]");//	add MSB to u, store in u.
		__asm("adcs t, t, #0x00");//	add carry to t, store in t.	

		temp1 = ((digit_t*)p503p1)[7];
		__asm("umull UV[0], UV[1], temp1, mc[10]");//	round 2: ii=17, jj=10;	multiply 8th dgt of p503p1 w/ 11th dgt of C.
		__asm("adds v, v, UV[0]");//	add LSB to v, store in v.
		__asm("adcs u, u, UV[1]");//	add MSB to u, store in u.
		__asm("adcs t, t, #0x00");//	add carry to t, store in t.
		
		__asm("adds mc[1], v, ma[17]");//	round 2: ii=17, jj=16;;	add v to 18th dgt of A, store in 2nd dgt of C.
		__asm("adcs v, u, #0x00");//	add carry to u, store in v.
		__asm("adcs u, t, #0x00");//	add carry to t, store in u.
		__asm("mov t, #0x00");//	set t to zero.


		temp1 = ((digit_t*)p503p1)[15];
		__asm("umull UV[0], UV[1], temp1, mc[3]");//	round 2: ii=18, jj=3;	multiply 16th dgt of p503p1 w/ 4th dgt of C.
		__asm("adds v, v, UV[0]");//	add LSB to v, store in v.
		__asm("adcs u, u, UV[1]");//	add MSB to u, store in u.
		__asm("adcs t, t, #0x00");//	add carry to t, store in t.
		
		temp1 = ((digit_t*)p503p1)[14];		
		__asm("umull UV[0], UV[1], temp1, mc[4]");//	round 2: ii=18, jj=4;	multiply 15th dgt of p503p1 w/ 5th dgt of C.
		__asm("adds v, v, UV[0]");//	add LSB to v, store in v.
		__asm("adcs u, u, UV[1]");//	add MSB to u, store in u.
		__asm("adcs t, t, #0x00");//	add carry to t, store in t.	

		temp1 = ((digit_t*)p503p1)[13];
		__asm("umull UV[0], UV[1], temp1, mc[5]");//	round 2: ii=18, jj=5;	multiply 14th dgt of p503p1 w/ 6th dgt of C.
		__asm("adds v, v, UV[0]");//	add LSB to v, store in v.
		__asm("adcs u, u, UV[1]");//	add MSB to u, store in u.
		__asm("adcs t, t, #0x00");//	add carry to t, store in t.
		
		temp1 = ((digit_t*)p503p1)[12];		
		__asm("umull UV[0], UV[1], temp1, mc[6]");//	round 2: ii=18, jj=6;	multiply 13th dgt of p503p1 w/ 7th dgt of C.
		__asm("adds v, v, UV[0]");//	add LSB to v, store in v.
		__asm("adcs u, u, UV[1]");//	add MSB to u, store in u.
		__asm("adcs t, t, #0x00");//	add carry to t, store in t.	

		temp1 = ((digit_t*)p503p1)[11];
		__asm("umull UV[0], UV[1], temp1, mc[7]");//	round 2: ii=18, jj=7;	multiply 12th dgt of p503p1 w/ 8th dgt of C.
		__asm("adds v, v, UV[0]");//	add LSB to v, store in v.
		__asm("adcs u, u, UV[1]");//	add MSB to u, store in u.
		__asm("adcs t, t, #0x00");//	add carry to t, store in t.
		
		temp1 = ((digit_t*)p503p1)[10];		
		__asm("umull UV[0], UV[1], temp1, mc[8]");//	round 2: ii=18, jj=8;	multiply 11th dgt of p503p1 w/ 9th dgt of C.
		__asm("adds v, v, UV[0]");//	add LSB to v, store in v.
		__asm("adcs u, u, UV[1]");//	add MSB to u, store in u.
		__asm("adcs t, t, #0x00");//	add carry to t, store in t.

		temp1 = ((digit_t*)p503p1)[9];
		__asm("umull UV[0], UV[1], temp1, mc[9]");//	round 2: ii=18, jj=9;	multiply 10th dgt of p503p1 w/ 10 dgt of C.
		__asm("adds v, v, UV[0]");//	add LSB to v, store in v.
		__asm("adcs u, u, UV[1]");//	add MSB to u, store in u.
		__asm("adcs t, t, #0x00");//	add carry to t, store in t.
		
		temp1 = ((digit_t*)p503p1)[8];		
		__asm("umull UV[0], UV[1], temp1, mc[10]");//	round 2: ii=18, jj=10;	multiply 9th dgt of p503p1 w/ 11th dgt of C.
		__asm("adds v, v, UV[0]");//	add LSB to v, store in v.
		__asm("adcs u, u, UV[1]");//	add MSB to u, store in u.
		__asm("adcs t, t, #0x00");//	add carry to t, store in t.

		temp1 = ((digit_t*)p503p1)[7];
		__asm("umull UV[0], UV[1], temp1, mc[11]");//	round 2: ii=18, jj=11;	multiply 8th dgt of p503p1 w/ 12th dgt of C.
		__asm("adds v, v, UV[0]");//	add LSB to v, store in v.
		__asm("adcs u, u, UV[1]");//	add MSB to u, store in u.
		__asm("adcs t, t, #0x00");//	add carry to t, store in t.
			
		__asm("adds mc[2], v, ma[18]");//	round 2: ii=18, jj=16;	add v to 19th dgt of A, store in 3rd dgt of C.
		__asm("adcs v, u, #0x00");//	add carry to u, store in v.
		__asm("adcs u, t, #0x00");//	add carry to t, store in u.
		__asm("mov t, #0x00");//	set t to zero.


		temp1 = ((digit_t*)p503p1)[15];
		__asm("umull UV[0], UV[1], temp1, mc[4]");//	round 2: ii=19, jj=4;	multiply 16th dgt of p503p1 w/ 5th dgt of C.
		__asm("adds v, v, UV[0]");//	add LSB to v, store in v.
		__asm("adcs u, u, UV[1]");//	add MSB to u, store in u.
		__asm("adcs t, t, #0x00");//	add carry to t, store in t.

		temp1 = ((digit_t*)p503p1)[14];
		__asm("umull UV[0], UV[1], temp1, mc[5]");//	round 2: ii=19, jj=5;	multiply 15th dgt of p503p1 w/ 6th dgt of C.
		__asm("adds v, v, UV[0]");//	add LSB to v, store in v.
		__asm("adcs u, u, UV[1]");//	add MSB to u, store in u.
		__asm("adcs t, t, #0x00");//	add carry to t, store in t.

		temp1 = ((digit_t*)p503p1)[13];
		__asm("umull UV[0], UV[1], temp1, mc[6]");//	round 2: ii=19, jj=6;	multiply 14th dgt of p503p1 w/ 7th dgt of C.
		__asm("adds v, v, UV[0]");//	add LSB to v, store in v.
		__asm("adcs u, u, UV[1]");//	add MSB to u, store in u.
		__asm("adcs t, t, #0x00");//	add carry to t, store in t.		

		temp1 = ((digit_t*)p503p1)[12];
		__asm("umull UV[0], UV[1], temp1, mc[7]");//	round 2: ii=19, jj=7;	multiply 13th dgt of p503p1 w/ 8th dgt of C.
		__asm("adds v, v, UV[0]");//	add LSB to v, store in v.
		__asm("adcs u, u, UV[1]");//	add MSB to u, store in u.
		__asm("adcs t, t, #0x00");//	add carry to t, store in t.		

		temp1 = ((digit_t*)p503p1)[11];
		__asm("umull UV[0], UV[1], temp1, mc[8]");//	round 2: ii=19, jj=8;	multiply 12th dgt of p503p1 w/ 9th dgt of C.
		__asm("adds v, v, UV[0]");//	add LSB to v, store in v.
		__asm("adcs u, u, UV[1]");//	add MSB to u, store in u.
		__asm("adcs t, t, #0x00");//	add carry to t, store in t.		
		
		temp1 = ((digit_t*)p503p1)[10];		
		__asm("umull UV[0], UV[1], temp1, mc[9]");//	round 2: ii=19, jj=9;	multiply 11th dgt of p503p1 w/ 10th dgt of C.
		__asm("adds v, v, UV[0]");//	add LSB to v, store in v.
		__asm("adcs u, u, UV[1]");//	add MSB to u, store in u.
		__asm("adcs t, t, #0x00");//	add carry to t, store in t.

		temp1 = ((digit_t*)p503p1)[9];
		__asm("umull UV[0], UV[1], temp1, mc[10]");//	round 2: ii=19, jj=10;	multiply 10th dgt of p503p1 w/ 11th dgt of C.
		__asm("adds v, v, UV[0]");//	add LSB to v, store in v.
		__asm("adcs u, u, UV[1]");//	add MSB to u, store in u.
		__asm("adcs t, t, #0x00");//	add carry to t, store in t.

		temp1 = ((digit_t*)p503p1)[8];
		__asm("umull UV[0], UV[1], temp1, mc[11]");//	round 2: ii=19, jj=11;	multiply 9th dgt of p503p1 w/ 12th dgt of C.
		__asm("adds v, v, UV[0]");//	add LSB to v, store in v.
		__asm("adcs u, u, UV[1]");//	add MSB to u, store in u.
		__asm("adcs t, t, #0x00");//	add carry to t, store in t.	

		temp1 = ((digit_t*)p503p1)[7];
		__asm("umull UV[0], UV[1], temp1, mc[12]");//	round 2: ii=19, jj=12;	multiply 8th dgt of p503p1 w/ 13th dgt of C.
		__asm("adds v, v, UV[0]");//	add LSB to v, store in v.
		__asm("adcs u, u, UV[1]");//	add MSB to u, store in u.
		__asm("adcs t, t, #0x00");//	add carry to t, store in t.
		
		__asm("adds mc[3], v, ma[19]");//	round 2: ii=19, jj=16;	add v to 20th dgt of A, store in 4th dgt of C.
		__asm("adcs v, u, #0x00");//	add carry to u, store in v.
		__asm("adcs u, t, #0x00");//	add carry to t, store in u.
		__asm("mov t, #0x00");//	set t to zero.


		temp1 = ((digit_t*)p503p1)[15];
		__asm("umull UV[0], UV[1], temp1, mc[5]");//	round 2: ii=20, jj=5;	multiply 16th dgt of p503p1 w/ 6th dgt of C.
		__asm("adds v, v, UV[0]");//	add LSB to v, store in v.
		__asm("adcs u, u, UV[1]");//	add MSB to u, store in u.
		__asm("adcs t, t, #0x00");//	add carry to t, store in t.
		
		temp1 = ((digit_t*)p503p1)[14];
		__asm("umull UV[0], UV[1], temp1, mc[6]");//	round 2: ii=20, jj=6;	multiply 15th dgt of p503p1 w/ 7th dgt of C.
		__asm("adds v, v, UV[0]");//	add LSB to v, store in v.
		__asm("adcs u, u, UV[1]");//	add MSB to u, store in u.
		__asm("adcs t, t, #0x00");//	add carry to t, store in t.		

		temp1 = ((digit_t*)p503p1)[13];
		__asm("umull UV[0], UV[1], temp1, mc[7]");//	round 2: ii=20, jj=7;	multiply 14th dgt of p503p1 w/ 8th dgt of C.
		__asm("adds v, v, UV[0]");//	add LSB to v, store in v.
		__asm("adcs u, u, UV[1]");//	add MSB to u, store in u.
		__asm("adcs t, t, #0x00");//	add carry to t, store in t.		

		temp1 = ((digit_t*)p503p1)[12];
		__asm("umull UV[0], UV[1], temp1, mc[8]");//	round 2: ii=20, jj=8;	multiply 13th dgt of p503p1 w/ 9th dgt of C.
		__asm("adds v, v, UV[0]");//	add LSB to v, store in v.
		__asm("adcs u, u, UV[1]");//	add MSB to u, store in u.
		__asm("adcs t, t, #0x00");//	add carry to t, store in t.		

		temp1 = ((digit_t*)p503p1)[11];
		__asm("umull UV[0], UV[1], temp1, mc[9]");//	round 2: ii=20, jj=9;	multiply 12th dgt of p503p1 w/ 10th dgt of C.
		__asm("adds v, v, UV[0]");//	add LSB to v, store in v.
		__asm("adcs u, u, UV[1]");//	add MSB to u, store in u.
		__asm("adcs t, t, #0x00");//	add carry to t, store in t.	

		temp1 = ((digit_t*)p503p1)[10];
		__asm("umull UV[0], UV[1], temp1, mc[10]");//	round 2: ii=20, jj=10;	multiply 11th dgt of p503p1 w/ 11th dgt of C.
		__asm("adds v, v, UV[0]");//	add LSB to v, store in v.
		__asm("adcs u, u, UV[1]");//	add MSB to u, store in u.
		__asm("adcs t, t, #0x00");//	add carry to t, store in t.
		
		temp1 = ((digit_t*)p503p1)[9];		
		__asm("umull UV[0], UV[1], temp1, mc[11]");//	round 2: ii=20, jj=11;	multiply 10th dgt of p503p1 w/ 12th dgt of C.
		__asm("adds v, v, UV[0]");//	add LSB to v, store in v.
		__asm("adcs u, u, UV[1]");//	add MSB to u, store in u.
		__asm("adcs t, t, #0x00");//	add carry to t, store in t.

		temp1 = ((digit_t*)p503p1)[8];
		__asm("umull UV[0], UV[1], temp1, mc[12]");//	round 2: ii=20, jj=12;	multiply 9th dgt of p503p1 w/ 13th dgt of C.
		__asm("adds v, v, UV[0]");//	add LSB to v, store in v.
		__asm("adcs u, u, UV[1]");//	add MSB to u, store in u.
		__asm("adcs t, t, #0x00");//	add carry to t, store in t.		
	
		temp1 = ((digit_t*)p503p1)[7];	
		__asm("umull UV[0], UV[1], temp1, mc[13]");//	round 2: ii=20, jj=13;	multiply 8th dgt of p503p1 w/ 14th dgt of C.
		__asm("adds v, v, UV[0]");//	add LSB to v, store in v.
		__asm("adcs u, u, UV[1]");//	add MSB to u, store in u.
		__asm("adcs t, t, #0x00");//	add carry to t, store in t.	

		__asm("adds mc[4], v, ma[20]");//	round 2: ii=20, jj=16;	add v to 21st dgt of A, store in 5th dgt of C.
		__asm("adcs v, u, #0x00");//	add carry to u, store in v.
		__asm("adcs u, t, #0x00");//	add carry to t, store in u.
		__asm("mov t, #0x00");//	set t to zero.
		

		temp1 = ((digit_t*)p503p1)[15];
		__asm("umull UV[0], UV[1], temp1, mc[6]");//	round 2: ii=21, jj=6;	multiply 16th dgt of p503p1 w/ 7th dgt of C.
		__asm("adds v, v, UV[0]");//	add LSB to v, store in v.
		__asm("adcs u, u, UV[1]");//	add MSB to u, store in u.
		__asm("adcs t, t, #0x00");//	add carry to t, store in t.	

		temp1 = ((digit_t*)p503p1)[14];
		__asm("umull UV[0], UV[1], temp1, mc[7]");//	round 2: ii=21, jj=7;	multiply 15th dgt of p503p1 w/ 8th dgt of C.
		__asm("adds v, v, UV[0]");//	add LSB to v, store in v.
		__asm("adcs u, u, UV[1]");//	add MSB to u, store in u.
		__asm("adcs t, t, #0x00");//	add carry to t, store in t.

		temp1 = ((digit_t*)p503p1)[13];
		__asm("umull UV[0], UV[1], temp1, mc[8]");//	round 2: ii=21, jj=8;	multiply 14th dgt of p503p1 w/ 9th dgt of C.
		__asm("adds v, v, UV[0]");//	add LSB to v, store in v.
		__asm("adcs u, u, UV[1]");//	add MSB to u, store in u.
		__asm("adcs t, t, #0x00");//	add carry to t, store in t.

		temp1 = ((digit_t*)p503p1)[12];
		__asm("umull UV[0], UV[1], temp1, mc[9]");//	round 2: ii=21, jj=9;	multiply 13th dgt of p503p1 w/ 10th dgt of C.
		__asm("adds v, v, UV[0]");//	add LSB to v, store in v.
		__asm("adcs u, u, UV[1]");//	add MSB to u, store in u.
		__asm("adcs t, t, #0x00");//	add carry to t, store in t.		

		temp1 = ((digit_t*)p503p1)[11];
		__asm("umull UV[0], UV[1], temp1, mc[10]");//	round 2: ii=21, jj=10;	multiply 12th dgt of p503p1 w/ 11th dgt of C.
		__asm("adds v, v, UV[0]");//	add LSB to v, store in v.
		__asm("adcs u, u, UV[1]");//	add MSB to u, store in u.
		__asm("adcs t, t, #0x00");//	add carry to t, store in t.		

		temp1 = ((digit_t*)p503p1)[10];
		__asm("umull UV[0], UV[1], temp1, mc[11]");//	round 2: ii=21, jj=11;	multiply 11th dgt of p503p1 w/ 12th dgt of C.
		__asm("adds v, v, UV[0]");//	add LSB to v, store in v.
		__asm("adcs u, u, UV[1]");//	add MSB to u, store in u.
		__asm("adcs t, t, #0x00");//	add carry to t, store in t.

		temp1 = ((digit_t*)p503p1)[9];
		__asm("umull UV[0], UV[1], temp1, mc[12]");//	round 2: ii=21, jj=12;	multiply 10th dgt of p503p1 w/ 13th dgt of C.
		__asm("adds v, v, UV[0]");//	add LSB to v, store in v.
		__asm("adcs u, u, UV[1]");//	add MSB to u, store in u.
		__asm("adcs t, t, #0x00");//	add carry to t, store in t.
		
		temp1 = ((digit_t*)p503p1)[8];		
		__asm("umull UV[0], UV[1], temp1, mc[13]");//	round 2: ii=21, jj=13;	multiply 9th dgt of p503p1 w/ 14th dgt of C.
		__asm("adds v, v, UV[0]");//	add LSB to v, store in v.
		__asm("adcs u, u, UV[1]");//	add MSB to u, store in u.
		__asm("adcs t, t, #0x00");//	add carry to t, store in t.

		temp1 = ((digit_t*)p503p1)[7];
		__asm("umull UV[0], UV[1], temp1, mc[14]");//	round 2: ii=21, jj=14;	multiply 8th dgt of p503p1 w/ 15th dgt of C.
		__asm("adds v, v, UV[0]");//	add LSB to v, store in v.
		__asm("adcs u, u, UV[1]");//	add MSB to u, store in u.
		__asm("adcs t, t, #0x00");//	add carry to t, store in t.
		
		__asm("adds mc[5], v, ma[21]");//	round 2: ii=21, jj=16;	add v to 22nd dgt of A, store in 6th dgt of C.
		__asm("adcs v, u, #0x00");//	add carry to u, store in v.
		__asm("adcs u, t, #0x00");//	add carry to t, store in u.
		__asm("mov t, #0x00");//	set t to zero.
		

		temp1 = ((digit_t*)p503p1)[15];
		__asm("umull UV[0], UV[1], temp1, mc[7]");//	round 2: ii=22, jj=7;	multiply 16th dgt of p503p1 w/ 8th dgt of C.
		__asm("adds v, v, UV[0]");//	add LSB to v, store in v.
		__asm("adcs u, u, UV[1]");//	add MSB to u, store in u.
		__asm("adcs t, t, #0x00");//	add carry to t, store in t.		

		temp1 = ((digit_t*)p503p1)[14];
		__asm("umull UV[0], UV[1], temp1, mc[8]");//	round 2: ii=22, jj=8;	multiply 15th dgt of p503p1 w/ 9th dgt of C.
		__asm("adds v, v, UV[0]");//	add LSB to v, store in v.
		__asm("adcs u, u, UV[1]");//	add MSB to u, store in u.
		__asm("adcs t, t, #0x00");//	add carry to t, store in t.	

		temp1 = ((digit_t*)p503p1)[13];
		__asm("umull UV[0], UV[1], temp1, mc[9]");//	round 2: ii=22, jj=9;	multiply 14th dgt of p503p1 w/ 10th dgt of C.
		__asm("adds v, v, UV[0]");//	add LSB to v, store in v.
		__asm("adcs u, u, UV[1]");//	add MSB to u, store in u.
		__asm("adcs t, t, #0x00");//	add carry to t, store in t.	

		temp1 = ((digit_t*)p503p1)[12];
		__asm("umull UV[0], UV[1], temp1, mc[10]");//	round 2: ii=22, jj=10;	multiply 13th dgt of p503p1 w/ 11th dgt of C.
		__asm("adds v, v, UV[0]");//	add LSB to v, store in v.
		__asm("adcs u, u, UV[1]");//	add MSB to u, store in u.
		__asm("adcs t, t, #0x00");//	add carry to t, store in t.		
		
		temp1 = ((digit_t*)p503p1)[11];		
		__asm("umull UV[0], UV[1], temp1, mc[11]");//	round 2: ii=22, jj=11;	multiply 12th dgt of p503p1 w/ 12th dgt of C.
		__asm("adds v, v, UV[0]");//	add LSB to v, store in v.
		__asm("adcs u, u, UV[1]");//	add MSB to u, store in u.
		__asm("adcs t, t, #0x00");//	add carry to t, store in t.		

		temp1 = ((digit_t*)p503p1)[10];
		__asm("umull UV[0], UV[1], temp1, mc[12]");//	round 2: ii=22, jj=12;	multiply 11th dgt of p503p1 w/ 13th dgt of C.
		__asm("adds v, v, UV[0]");//	add LSB to v, store in v.
		__asm("adcs u, u, UV[1]");//	add MSB to u, store in u.
		__asm("adcs t, t, #0x00");//	add carry to t, store in t.		
		
		temp1 = ((digit_t*)p503p1)[9];		
		__asm("umull UV[0], UV[1], temp1, mc[13]");//	round 2: ii=22, jj=13;	multiply 10th dgt of p503p1 w/ 14th dgt of C.
		__asm("adds v, v, UV[0]");//	add LSB to v, store in v.
		__asm("adcs u, u, UV[1]");//	add MSB to u, store in u.
		__asm("adcs t, t, #0x00");//	add carry to t, store in t.	

		temp1 = ((digit_t*)p503p1)[8];
		__asm("umull UV[0], UV[1], temp1, mc[14]");//	round 2: ii=22, jj=14;	multiply 9th dgt of p503p1 w/ 15th dgt of C.
		__asm("adds v, v, UV[0]");//	add LSB to v, store in v.
		__asm("adcs u, u, UV[1]");//	add MSB to u, store in u.
		__asm("adcs t, t, #0x00");//	add carry to t, store in t.	

		temp1 = ((digit_t*)p503p1)[7];
		__asm("umull UV[0], UV[1], temp1, mc[15]");//	round 2: ii=22, jj=15;	multiply 8th dgt of p503p1 w/ 16th dgt of C.
		__asm("adds v, v, UV[0]");//	add LSB to v, store in v.
		__asm("adcs u, u, UV[1]");//	add MSB to u, store in u.
		__asm("adcs t, t, #0x00");//	add carry to t, store in t.	

		__asm("adds mc[6], v, ma[22]");//	round 2: ii=22, jj=16;	add v to 23rd dgt of A, store in 7th dgt of C.
		__asm("adcs v, u, #0x00");//	add carry to u, store in v.
		__asm("adcs u, t, #0x00");//	add carry to t, store in u.
		__asm("mov t, #0x00");//	set t to zero.
		
		
		temp1 = ((digit_t*)p503p1)[15];		
		__asm("umull UV[0], UV[1], temp1, mc[8]");//	round 2: ii=23, jj=8;	multiply 16th dgt of p503p1 w/ 9th dgt of C.
		__asm("adds v, v, UV[0]");//	add LSB to v, store in v.
		__asm("adcs u, u, UV[1]");//	add MSB to u, store in u.
		__asm("adcs t, t, #0x00");//	add carry to t, store in t.

		temp1 = ((digit_t*)p503p1)[14];
		__asm("umull UV[0], UV[1], temp1, mc[9]");//	round 2: ii=23, jj=9;	multiply 15th dgt of p503p1 w/ 10 dgt of C.
		__asm("adds v, v, UV[0]");//	add LSB to v, store in v.
		__asm("adcs u, u, UV[1]");//	add MSB to u, store in u.
		__asm("adcs t, t, #0x00");//	add carry to t, store in t.	

		temp1 = ((digit_t*)p503p1)[13];
		__asm("umull UV[0], UV[1], temp1, mc[10]");//	round 2: ii=23, jj=10;	multiply 14th dgt of p503p1 w/ 11th dgt of C.
		__asm("adds v, v, UV[0]");//	add LSB to v, store in v.
		__asm("adcs u, u, UV[1]");//	add MSB to u, store in u.
		__asm("adcs t, t, #0x00");//	add carry to t, store in t.		

		temp1 = ((digit_t*)p503p1)[12];
		__asm("umull UV[0], UV[1], temp1, mc[11]");//	round 2: ii=23, jj=11;	multiply 13th dgt of p503p1 w/ 12th dgt of C.
		__asm("adds v, v, UV[0]");//	add LSB to v, store in v.
		__asm("adcs u, u, UV[1]");//	add MSB to u, store in u.
		__asm("adcs t, t, #0x00");//	add carry to t, store in t.	

		temp1 = ((digit_t*)p503p1)[11];
		__asm("umull UV[0], UV[1], temp1, mc[12]");//	round 2: ii=23, jj=12;	multiply 12th dgt of p503p1 w/ 13th dgt of C.
		__asm("adds v, v, UV[0]");//	add LSB to v, store in v.
		__asm("adcs u, u, UV[1]");//	add MSB to u, store in u.
		__asm("adcs t, t, #0x00");//	add carry to t, store in t.			

		temp1 = ((digit_t*)p503p1)[10];
		__asm("umull UV[0], UV[1], temp1, mc[13]");//	round 2: ii=23, jj=13;	multiply 11th dgt of p503p1 w/ 14th dgt of C.
		__asm("adds v, v, UV[0]");//	add LSB to v, store in v.
		__asm("adcs u, u, UV[1]");//	add MSB to u, store in u.
		__asm("adcs t, t, #0x00");//	add carry to t, store in t.	

		temp1 = ((digit_t*)p503p1)[9];
		__asm("umull UV[0], UV[1], temp1, mc[14]");//	round 2: ii=23, jj=14;	multiply 10th dgt of p503p1 w/ 15 dgt of C.
		__asm("adds v, v, UV[0]");//	add LSB to v, store in v.
		__asm("adcs u, u, UV[1]");//	add MSB to u, store in u.
		__asm("adcs t, t, #0x00");//	add carry to t, store in t.	

		temp1 = ((digit_t*)p503p1)[8];
		__asm("umull UV[0], UV[1], temp1, mc[15]");//	round 2: ii=23, jj=15;	multiply 9th dgt of p503p1 w/ 16th dgt of C.
		__asm("adds v, v, UV[0]");//	add LSB to v, store in v.
		__asm("adcs u, u, UV[1]");//	add MSB to u, store in u.
		__asm("adcs t, t, #0x00");//	add carry to t, store in t.

		__asm("adds mc[7], v, ma[23]");//	round 2: ii=23, jj=16;	add v to 24th dgt of A, store in 8th dgt of C.
		__asm("adcs v, u, #0x00");//	add carry to u, store in v.
		__asm("adcs u, t, #0x00");//	add carry to t, store in u.
		__asm("mov t, #0x00");//	set t to zero.
		

		temp1 = ((digit_t*)p503p1)[15];
		__asm("umull UV[0], UV[1], temp1, mc[9]");//	round 2: ii=24, jj=9;	multiply 16th dgt of p503p1 w/ 10th dgt of C.
		__asm("adds v, v, UV[0]");//	add LSB to v, store in v.
		__asm("adcs u, u, UV[1]");//	add MSB to u, store in u.
		__asm("adcs t, t, #0x00");//	add carry to t, store in t.		

		temp1 = ((digit_t*)p503p1)[14];
		__asm("umull UV[0], UV[1], temp1, mc[10]");//	round 2: ii=24, jj=10;	multiply 15th dgt of p503p1 w/ 11th dgt of C.
		__asm("adds v, v, UV[0]");//	add LSB to v, store in v.
		__asm("adcs u, u, UV[1]");//	add MSB to u, store in u.
		__asm("adcs t, t, #0x00");//	add carry to t, store in t.			

		temp1 = ((digit_t*)p503p1)[13];
		__asm("umull UV[0], UV[1], temp1, mc[11]");//	round 2: ii=24, jj=11;	multiply 14th dgt of p503p1 w/ 12th dgt of C.
		__asm("adds v, v, UV[0]");//	add LSB to v, store in v.
		__asm("adcs u, u, UV[1]");//	add MSB to u, store in u.
		__asm("adcs t, t, #0x00");//	add carry to t, store in t.	

		temp1 = ((digit_t*)p503p1)[12];
		__asm("umull UV[0], UV[1], temp1, mc[12]");//	round 2: ii=24, jj=12;	multiply 13th dgt of p503p1 w/ 13th dgt of C.
		__asm("adds v, v, UV[0]");//	add LSB to v, store in v.
		__asm("adcs u, u, UV[1]");//	add MSB to u, store in u.
		__asm("adcs t, t, #0x00");//	add carry to t, store in t.	

		temp1 = ((digit_t*)p503p1)[11];
		__asm("umull UV[0], UV[1], temp1, mc[13]");//	round 2: ii=24, jj=13;	multiply 12th dgt of p503p1 w/ 14th dgt of C.
		__asm("adds v, v, UV[0]");//	add LSB to v, store in v.
		__asm("adcs u, u, UV[1]");//	add MSB to u, store in u.
		__asm("adcs t, t, #0x00");//	add carry to t, store in t.		

		temp1 = ((digit_t*)p503p1)[10];
		__asm("umull UV[0], UV[1], temp1, mc[14]");//	round 2: ii=24, jj=14;	multiply 11th dgt of p503p1 w/ 15th dgt of C.
		__asm("adds v, v, UV[0]");//	add LSB to v, store in v.
		__asm("adcs u, u, UV[1]");//	add MSB to u, store in u.
		__asm("adcs t, t, #0x00");//	add carry to t, store in t.	

		temp1 = ((digit_t*)p503p1)[9];
		__asm("umull UV[0], UV[1], temp1, mc[15]");//	round 2: ii=24, jj=15;	multiply 10th dgt of p503p1 w/ 16th dgt of C.
		__asm("adds v, v, UV[0]");//	add LSB to v, store in v.
		__asm("adcs u, u, UV[1]");//	add MSB to u, store in u.
		__asm("adcs t, t, #0x00");//	add carry to t, store in t.

		__asm("adds mc[8], v, ma[24]");//	round 2: ii=24, jj=16;	add v to 25th dgt of A, store in 9th dgt of C.
		__asm("adcs v, u, #0x00");//	add carry to u, store in v.
		__asm("adcs u, t, #0x00");//	add carry to t, store in u.
		__asm("mov t, #0x00");//	set t to zero.
		
		
		temp1 = ((digit_t*)p503p1)[15];		
		__asm("umull UV[0], UV[1], temp1, mc[10]");//	round 2: ii=25, jj=10;	multiply 16th dgt of p503p1 w/ 11th dgt of C.
		__asm("adds v, v, UV[0]");//	add LSB to v, store in v.
		__asm("adcs u, u, UV[1]");//	add MSB to u, store in u.
		__asm("adcs t, t, #0x00");//	add carry to t, store in t.		
		
		temp1 = ((digit_t*)p503p1)[14];		
		__asm("umull UV[0], UV[1], temp1, mc[11]");//	round 2: ii=25, jj=11;	multiply 15th dgt of p503p1 w/ 12th dgt of C.
		__asm("adds v, v, UV[0]");//	add LSB to v, store in v.
		__asm("adcs u, u, UV[1]");//	add MSB to u, store in u.
		__asm("adcs t, t, #0x00");//	add carry to t, store in t.		

		temp1 = ((digit_t*)p503p1)[13];
		__asm("umull UV[0], UV[1], temp1, mc[12]");//	round 2: ii=25, jj=12;	multiply 14th dgt of p503p1 w/ 13th dgt of C.
		__asm("adds v, v, UV[0]");//	add LSB to v, store in v.
		__asm("adcs u, u, UV[1]");//	add MSB to u, store in u.
		__asm("adcs t, t, #0x00");//	add carry to t, store in t.	

		temp1 = ((digit_t*)p503p1)[12];
		__asm("umull UV[0], UV[1], temp1, mc[13]");//	round 2: ii=25, jj=13;	multiply 13th dgt of p503p1 w/ 14th dgt of C.
		__asm("adds v, v, UV[0]");//	add LSB to v, store in v.
		__asm("adcs u, u, UV[1]");//	add MSB to u, store in u.
		__asm("adcs t, t, #0x00");//	add carry to t, store in t.
		
		temp1 = ((digit_t*)p503p1)[11];		
		__asm("umull UV[0], UV[1], temp1, mc[14]");//	round 2: ii=25, jj=14;	multiply 12th dgt of p503p1 w/ 15th dgt of C.
		__asm("adds v, v, UV[0]");//	add LSB to v, store in v.
		__asm("adcs u, u, UV[1]");//	add MSB to u, store in u.
		__asm("adcs t, t, #0x00");//	add carry to t, store in t.	

		temp1 = ((digit_t*)p503p1)[10];
		__asm("umull UV[0], UV[1], temp1, mc[15]");//	round 2: ii=25, jj=15;	multiply 11th dgt of p503p1 w/ 16th dgt of C.
		__asm("adds v, v, UV[0]");//	add LSB to v, store in v.
		__asm("adcs u, u, UV[1]");//	add MSB to u, store in u.
		__asm("adcs t, t, #0x00");//	add carry to t, store in t.
			
		__asm("adds mc[9], v, ma[25]");//	round 2: ii=25, jj=16;	add v to 26th dgt of A, store in 10th dgt of C.
		__asm("adcs v, u, #0x00");//	add carry to u, store in v.
		__asm("adcs u, t, #0x00");//	add carry to t, store in u.
		__asm("mov t, #0x00");//	set t to zero.		
		
		
		temp1 = ((digit_t*)p503p1)[15];		
		__asm("umull UV[0], UV[1], temp1, mc[11]");//	round 2: ii=26, jj=11;	multiply 16th dgt of p503p1 w/ 12th dgt of C.
		__asm("adds v, v, UV[0]");//	add LSB to v, store in v.
		__asm("adcs u, u, UV[1]");//	add MSB to u, store in u.
		__asm("adcs t, t, #0x00");//	add carry to t, store in t.	

		temp1 = ((digit_t*)p503p1)[14];
		__asm("umull UV[0], UV[1], temp1, mc[12]");//	round 2: ii=26, jj=12;	multiply 15th dgt of p503p1 w/ 13th dgt of C.
		__asm("adds v, v, UV[0]");//	add LSB to v, store in v.
		__asm("adcs u, u, UV[1]");//	add MSB to u, store in u.
		__asm("adcs t, t, #0x00");//	add carry to t, store in t.	

		temp1 = ((digit_t*)p503p1)[13];
		__asm("umull UV[0], UV[1], temp1, mc[13]");//	round 2: ii=26, jj=13;	multiply 14th dgt of p503p1 w/ 14th dgt of C.
		__asm("adds v, v, UV[0]");//	add LSB to v, store in v.
		__asm("adcs u, u, UV[1]");//	add MSB to u, store in u.
		__asm("adcs t, t, #0x00");//	add carry to t, store in t.			
		
		temp1 = ((digit_t*)p503p1)[12];		
		__asm("umull UV[0], UV[1], temp1, mc[14]");//	round 2: ii=26, jj=14;	multiply 13th dgt of p503p1 w/ 15th dgt of C.
		__asm("adds v, v, UV[0]");//	add LSB to v, store in v.
		__asm("adcs u, u, UV[1]");//	add MSB to u, store in u.
		__asm("adcs t, t, #0x00");//	add carry to t, store in t.			
		
		temp1 = ((digit_t*)p503p1)[11];		
		__asm("umull UV[0], UV[1], temp1, mc[15]");//	round 2: ii=26, jj=15;	multiply 12th dgt of p503p1 w/ 16th dgt of C.
		__asm("adds v, v, UV[0]");//	add LSB to v, store in v.
		__asm("adcs u, u, UV[1]");//	add MSB to u, store in u.
		__asm("adcs t, t, #0x00");//	add carry to t, store in t.		

		__asm("adds mc[10], v, ma[26]");//	round 2: ii=26, jj=16;	add v to 27th dgt of A, store in 11th dgt of C.
		__asm("adcs v, u, #0x00");//	add carry to u, store in v.
		__asm("adcs u, t, #0x00");//	add carry to t, store in u.
		__asm("mov t, #0x00");//	set t to zero.
		

		temp1 = ((digit_t*)p503p1)[15];
		__asm("umull UV[0], UV[1], temp1, mc[12]");//	round 2: ii=27, jj=12;	multiply 16th dgt of p503p1 w/ 13th dgt of C.
		__asm("adds v, v, UV[0]");//	add LSB to v, store in v.
		__asm("adcs u, u, UV[1]");//	add MSB to u, store in u.
		__asm("adcs t, t, #0x00");//	add carry to t, store in t.	

		temp1 = ((digit_t*)p503p1)[14];
		__asm("umull UV[0], UV[1], temp1, mc[13]");//	round 2: ii=27, jj=13;	multiply 15th dgt of p503p1 w/ 14th dgt of C.
		__asm("adds v, v, UV[0]");//	add LSB to v, store in v.
		__asm("adcs u, u, UV[1]");//	add MSB to u, store in u.
		__asm("adcs t, t, #0x00");//	add carry to t, store in t.	

		temp1 = ((digit_t*)p503p1)[13];
		__asm("umull UV[0], UV[1], temp1, mc[14]");//	round 2: ii=27, jj=14;	multiply 14th dgt of p503p1 w/ 15th dgt of C.
		__asm("adds v, v, UV[0]");//	add LSB to v, store in v.
		__asm("adcs u, u, UV[1]");//	add MSB to u, store in u.
		__asm("adcs t, t, #0x00");//	add carry to t, store in t.	
		
		temp1 = ((digit_t*)p503p1)[12];		
		__asm("umull UV[0], UV[1], temp1, mc[15]");//	round 2: ii=27, jj=15;	multiply 13th dgt of p503p1 w/ 16th dgt of C.
		__asm("adds v, v, UV[0]");//	add LSB to v, store in v.
		__asm("adcs u, u, UV[1]");//	add MSB to u, store in u.
		__asm("adcs t, t, #0x00");//	add carry to t, store in t.		

		__asm("adds mc[11], v, ma[27]");//	round 2: ii=27, jj=16;	add v to 28th dgt of A, store in 12th dgt of C.
		__asm("adcs v, u, #0x00");//	add carry to u, store in v.
		__asm("adcs u, t, #0x00");//	add carry to t, store in u.
		__asm("mov t, #0x00");//	set t to zero.
		
		
		temp1 = ((digit_t*)p503p1)[15];		
		__asm("umull UV[0], UV[1], temp1, mc[13]");//	round 2: ii=28, jj=13;	multiply 16th dgt of p503p1 w/ 14th dgt of C.
		__asm("adds v, v, UV[0]");//	add LSB to v, store in v.
		__asm("adcs u, u, UV[1]");//	add MSB to u, store in u.
		__asm("adcs t, t, #0x00");//	add carry to t, store in t.	

		temp1 = ((digit_t*)p503p1)[14];
		__asm("umull UV[0], UV[1], temp1, mc[14]");//	round 2: ii=28, jj=14;	multiply 15th dgt of p503p1 w/ 15th dgt of C.
		__asm("adds v, v, UV[0]");//	add LSB to v, store in v.
		__asm("adcs u, u, UV[1]");//	add MSB to u, store in u.
		__asm("adcs t, t, #0x00");//	add carry to t, store in t.
		
		temp1 = ((digit_t*)p503p1)[13];		
		__asm("umull UV[0], UV[1], temp1, mc[15]");//	round 2: ii=28, jj=15;	multiply 14th dgt of p503p1 w/ 16th dgt of C.
		__asm("adds v, v, UV[0]");//	add LSB to v, store in v.
		__asm("adcs u, u, UV[1]");//	add MSB to u, store in u.
		__asm("adcs t, t, #0x00");//	add carry to t, store in t.		
		
		__asm("adds mc[12], v, ma[28]");//	round 2: ii=28, jj=16;	add v to 29th dgt of A, store in 13th dgt of C.
		__asm("adcs v, u, #0x00");//	add carry to u, store in v.
		__asm("adcs u, t, #0x00");//	add carry to t, store in u.
		__asm("mov t, #0x00");//	set t to zero.


		temp1 = ((digit_t*)p503p1)[15];
		__asm("umull UV[0], UV[1], temp1, mc[14]");//	round 2: ii=29, jj=14;	multiply 16th dgt of p503p1 w/ 15th dgt of C.
		__asm("adds v, v, UV[0]");//	add LSB to v, store in v.
		__asm("adcs u, u, UV[1]");//	add MSB to u, store in u.
		__asm("adcs t, t, #0x00");//	add carry to t, store in t.		
		
		temp1 = ((digit_t*)p503p1)[14];		
		__asm("umull UV[0], UV[1], temp1, mc[15]");//	round 2: ii=29, jj=15;	multiply 15th dgt of p503p1 w/ 16th dgt of C.
		__asm("adds v, v, UV[0]");//	add LSB to v, store in v.
		__asm("adcs u, u, UV[1]");//	add MSB to u, store in u.
		__asm("adcs t, t, #0x00");//	add carry to t, store in t.

		__asm("adds mc[13], v, ma[29]");//	round 2: ii=29, jj=16;	add v to 30th dgt of A, store in 14th dgt of C.
		__asm("adcs v, u, #0x00");//	add carry to u, store in v.
		__asm("adcs u, t, #0x00");//	add carry to t, store in u.
		__asm("mov t, #0x00");//	set t to zero.
		
		
		temp1 = ((digit_t*)p503p1)[15];		
		__asm("umull UV[0], UV[1], temp1, mc[15]");//	round 2: ii=30, jj=15;	multiply 16th dgt of p503p1 w/ 16th dgt of C.
		__asm("adds v, v, UV[0]");//	add LSB to v, store in v.
		__asm("adcs u, u, UV[1]");//	add MSB to u, store in u.
		__asm("adcs t, t, #0x00");//	add carry to t, store in t.	

		__asm("adds mc[14], v, ma[30]");//	round 2: ii=30, jj=16;	add v to 31st dgt of A, store in 15th dgt of C.
		__asm("adcs v, u, #0x00");//	add carry to u, store in v.
		__asm("adcs u, t, #0x00");//	add carry to t, store in u.
		__asm("mov t, #0x00");//	set t to zero.
		
		__asm("adds mc[15], v, ma[31]");
		
//	
//    for (i = 0; i < NWORDS_FIELD; i++) {
//        for (j = 0; j < i; j++) {
//            if (j < (i-p503_ZERO_WORDS+1)) { 
//							temp1 = ((digit_t*)p503p1)[i-j];
//               // MUL(mc[j], ((digit_t*)p503p1)[i-j], UV+1, UV[0]);
//								__asm("umull UV[0],UV[1],mc[j],temp1");
//                //ADDC(0, UV[0], v, carry, v); 
//                //ADDC(carry, UV[1], u, carry, u); 
//                //t += carry; 
//							__asm("adds v,v,UV[0]");
//							__asm("adcs u,u,UV[1]");
//							__asm("adcs t,t,#0x00");
//            }
//        }
//        //ADDC(0, v, ma[i], carry, v); 
//        //ADDC(carry, u, 0, carry, u); 				
//       //t += carry; 
//				__asm("adds v,v,ma[i]");
//				__asm("adcs u,u,#0x00");
//				__asm("adcs t,t,#0x00");
//        mc[i] = v;
//        v = u;
//        u = t;
//        t = 0;
//    }    

//    for (i = NWORDS_FIELD; i < 2*NWORDS_FIELD-1; i++) {
//        if (count > 0) {
//            count -= 1;
//        }
//        for (j = i-NWORDS_FIELD+1; j < NWORDS_FIELD; j++) {
//            if (j < (NWORDS_FIELD-count)) { 
//							temp1 = ((digit_t*)p503p1)[i-j];
//              //  MUL(mc[j], ((digit_t*)p503p1)[i-j], UV+1, UV[0]);
//							__asm("umull UV[0],UV[1],mc[j],temp1");
//                //ADDC(0, UV[0], v, carry, v); 
//                //ADDC(carry, UV[1], u, carry, u); 
//                //t += carry;
//							__asm("adds v,v,UV[0]");
//							__asm("adcs u,u,UV[1]");
//							__asm("adcs t,t,#0x00");
//            }
//        }
//        //ADDC(0, v, ma[i], carry, v); 
//        //ADDC(carry, u, 0, carry, u); 
//        //t += carry; 
//				__asm("adds v,v,ma[i]");
//				__asm("adcs u,u,#0x00");
//				__asm("adcs t,t,#0x00");
//        mc[i-NWORDS_FIELD] = v;
//        v = u;
//        u = t;
//        t = 0;
//    }
//    //ADDC(0, v, ma[2*NWORDS_FIELD-1], carry, v); 
//		__asm("adds v,v,ma[31]");		// ma[2*NWORDS_FIELD-1]
//    //mc[NWORDS_FIELD-1] = v;
//		mc[15] = v;
}
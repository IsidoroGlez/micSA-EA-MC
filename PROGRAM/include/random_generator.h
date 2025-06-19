#include <stdint.h>

//A parte de lo que valga NBR, en algunos sitios queremos de verdad 64 bits:

#define TWO64 ((double) 18446744073709551616.)   //2^64 
#define TWO64MINUSONE 18446744073709551615ULL //2^64 - 1
#define TWO63 9223372036854775808ULL  //2^63 

//Y tambien 24 bits para el generador de Luscher:
#define TWO24       16777216
#define TWO24M1     16777215

#if(NBR==32)
#define FNORM (2.32830636E-10F)// max float such that RAND_MAX*NORMF<1
#define TWOBRMINUS1   4294967295
#define TWOBR         ((double) 4294967296.)
#define TWOBRM1       2147483648U
#define FNORM (2.32830636E-10F)// max double such that RAND_MAX*NORMF<1
typedef uint32_t randint;
#else
#define TWOBRMINUS1   18446744073709551615ULL //2^64 - 1
#define TWOBR         ((double) 18446744073709551616.)   //2^64 
#define TWOBRM1       9223372036854775808ULL  //2^63 
#define FNORM (5.4210108624275218e-20)// max double such that RAND_MAX*NORMF<1
typedef uint64_t randint;
#endif

typedef struct   //Data_structure describing a good 64 bits generator
{
  uint64_t ira[256];
  uint64_t PR;
  uint64_t zseed;
  uint64_t final;
  unsigned char ip,ip1,ip2,ip3;
  unsigned char bobo[4];
} s_aleatorio_HQ_64bits;


typedef struct   //Data_structure describing a good 64 bits generator
{
  uint64_t s[4];
  uint64_t t;
  uint64_t final;
} s_xoshiro256pp;


#define _rotate_left(x,k) ((x << k) | (x >> (64 - k)))


#define _actualiza_xoshiro256pp(u){			\
    u.final=_rotate_left((u.s[0]+u.s[3]),23)+u.s[0];    \
    u.t=u.s[1]<<17;					\
    u.s[2]^=u.s[0];					\
    u.s[3]^=u.s[1];					\
    u.s[1]^=u.s[2];					\
    u.s[0]^=u.s[3];					\
    u.s[2]^=u.t;					\
    u.s[3]=_rotate_left(u.s[3], 45);			\
}


#define _actualiza_aleatorio_HQ_escalar(u) {\
        u.ira[u.ip]=u.ira[u.ip1]+u.ira[u.ip2];  \
        u.PR=u.ira[u.ip]^u.ira[u.ip3];          \
        u.zseed=u.zseed*3202034522624059733LLU+1;       \
        u.final=u.zseed+u.PR;                           \
        u.ip++;u.ip1++;u.ip2++;u.ip3++;}


void Init_Rand_HQ_64bits(s_aleatorio_HQ_64bits *, uint64_t);

#ifdef JANUS
#ifdef MAIN
randint zseed;
unsigned char ip, ip1, ip2, ip3;
randint ira[256];
#else
extern randint zseed;
extern unsigned char ip, ip1, ip2, ip3;
extern randint ira[256];
#endif

#define CGRANDOM  ( zseed=zseed*3202034522624059733LLU+1)
#define PRRANDOM  ( (ira[ip++]=ira[ip1++]+ira[ip2++]) ^ira[ip3++] )
#define HIRANDOM  ((randint) (PRRANDOM + CGRANDOM)) //Addition of both
void ini_random(unsigned long long);
#endif

uint64_t comprueba_semilla(uint64_t);
uint64_t lee_urandom(void);
void Init_luxury(uint64_t);
uint32_t luxury(void);//Generador de Luescher
uint64_t Luescher_PRNG(void);//8 bytes de 8 llamadas consecutivas a luxury.
void Inicia_generadores_CPU(s_aleatorio_HQ_64bits *, s_xoshiro256pp *,
			    uint32_t *, uint64_t);
void generate_seeds_from_one(uint64_t *, uint64_t *, int);

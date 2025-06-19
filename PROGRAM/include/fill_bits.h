/* MAIN: (for convenience, the placement of the counter bits is done in a separate macro)                                                                                                                   
The 40 time bits are distributed in (15+parity bit) in .x and 25 in .y                                                                                                                                     
In .x we have (at most) 2 internal counter bits and (15+1) external counter bits. That is 18 bits.                                                                                                             
There is space for 14 bits of entropy.                                                                                                                                                                            
In .y we have 25 bits of external counter, we have space for 7 bits of entropy.                                                                                                                              
In .z we have 18 bits reserved for the thread index, tid. There is space for 14 bits of entropy.                                                                                                                       
In .w we have 14 bits reserved for the jobID index. We have space left for 18 bits of entropy */

#define NUMBITS_TOTALTHREADS_IN_METROPOLIS  19
#define NUMBITS_PER_JOBID  14
#define MAX_BITS_RUN_LENGTH_IN_METROPOLIS_STEPS 40
#define NBITS_INTERNAL_CALL  2   //2^NBITS_INTERNAL_CALL must be >= number of calls to Salmdons in the kernel

//Checks y defines derivados

#if(SNBETASNR>(1<<NUMBITS_TOTALTHREADS_IN_METROPOLIS))
#error "Recompilar aumentando NUMBITS_TOTALTHREADS_IN_METROPOLIS"
#endif


#define TOTAL_BITS_COUNTER (NUMBITS_PER_JOBID+NBITS_INTERNAL_CALL+MAX_BITS_RUN_LENGTH_IN_METROPOLIS_STEPS+NUMBITS_TOTALTHREADS_IN_METROPOLIS+1)
#define NUM_ENTROPY_BITS (128 - TOTAL_BITS_COUNTER)
#if(NUM_ENTROPY_BITS>64)
#error "We need to generate more than 64 bits of entropy per call!"
#endif
#if(NUM_ENTROPY_BITS<0)
#error "128-bits counter too short!"
#endif

#if(NBITS_INTERNAL_CALL>16)
#error "Not enough space for internal calls in .x word!"
#endif
#if(MAX_BITS_RUN_LENGTH_IN_METROPOLIS_STEPS<15)
#error "We need at least 16 bits in time duration"
#endif
#if(MAX_BITS_RUN_LENGTH_IN_METROPOLIS_STEPS>47)
#error "This scheme is good for total time <= 2^{47}-1"
#endif

#define NBITS_ENTROPY_MASK_X (16-NBITS_INTERNAL_CALL)
#define ENTROPY_MASK_X ((1ULL<<NBITS_ENTROPY_MASK_X)-1ULL)
#define NBITS_ENTROPY_MASK_Y (47 - MAX_BITS_RUN_LENGTH_IN_METROPOLIS_STEPS) 
#define ENTROPY_MASK_Y ((1ULL<<NBITS_ENTROPY_MASK_Y)-1ULL)
#define NBITS_ENTROPY_MASK_Z (32 - NUMBITS_TOTALTHREADS_IN_METROPOLIS) 
#define ENTROPY_MASK_Z ((1ULL<<NBITS_ENTROPY_MASK_Z)-1ULL)

#define _fill_bits {                                                    \
  for(ir=0;ir<NR;ir++){							\
    _actualiza_xoshiro256pp(random_xoshiro256pp[ir]);			\
    _actualiza_aleatorio_HQ_escalar(random_PRC[ir]);			\
    entropy=random_PRC[ir].final+random_xoshiro256pp[ir].final;		\
    entropy>>=64-NUM_ENTROPY_BITS;                                      \
    temporal=iter>>15;                                                  \
    s_time_and_entropy.vec[ir].y=(uint32_t) temporal;			\
    temporal=(iter&32767ULL)<<(NBITS_INTERNAL_CALL+1);			\
    temporal|=(parity&1ULL)<<NBITS_INTERNAL_CALL;			\
    s_time_and_entropy.vec[ir].x=(uint32_t) temporal;			\
    temporal=(entropy&ENTROPY_MASK_X)<<(32-NBITS_ENTROPY_MASK_X);       \
    s_time_and_entropy.vec[ir].x|=(uint32_t) temporal;			\
    entropy>>=NBITS_ENTROPY_MASK_X;                                     \
    temporal=(entropy&ENTROPY_MASK_Y)<<(32-NBITS_ENTROPY_MASK_Y);       \
    s_time_and_entropy.vec[ir].y|=(uint32_t) temporal;			\
    entropy>>=NBITS_ENTROPY_MASK_Y;                                     \
    temporal=(entropy&ENTROPY_MASK_Z)<<(32-NBITS_ENTROPY_MASK_Z);       \
    s_time_and_entropy.vec[ir].z=(uint32_t) temporal;			\
    entropy>>=NBITS_ENTROPY_MASK_Z;                                     \
    temporal=entropy<<NUMBITS_PER_JOBID;                                \
    s_time_and_entropy.vec[ir].w=((uint32_t) temporal)|mi_jobID;}}

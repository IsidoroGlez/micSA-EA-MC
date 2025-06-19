#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <stdarg.h>
#include <sys/types.h>
#include <quadmath.h>

#include "random_generator.h"

#define ISI

// DEFINITIONS

#define MAXSTR 1024
#define NR 16
#define WARPSIZE 32
#define CUDATHREADS 512
#define MAXNSAMPLES 2

// CHECK DEFINES

#ifndef NUMBITSPREBUSQUEDAS
#error "NUMBITSPREBUSQUEDAS must be defined"
#endif

#define NPREBUSQUEDAS (1<<NUMBITSPREBUSQUEDAS)

#ifndef NBETAS
#error "NBETAS must be defined"
#endif

#ifndef L
#error "L must be defined"
#endif

#if(((L>>3)<<3)!=L)
#error "L must be a multiple of 8"
#endif

// PHYSICAL SYSTEM VARIABLES
#define MAXNWALKERS 12

#define DIM 3
#define DEGREE (DIM<<1)
#define Lx L
#define Ly L
#define Lz L
#define V (Lx*Ly*Lz)
#define S (Ly*Lx)

#define NRNBETAS (NR*NBETAS)

//MSC SYSTEM VARIABLES
#define MSC_L (L>>2)
#define MSC_Lz (Lz>>2)
#define MSC_S (MSC_L*MSC_L)
#define MSC_V (MSC_S*MSC_Lz)
#define HALF_MSC_V (MSC_V>>1)
#define HALF_MSC_VNR (HALF_MSC_V*NR)
#define HALF_MSC_VNRNBETAS (HALF_MSC_VNR*NBETAS)
#define MSC_VDEGREE (MSC_V*DEGREE)

enum{MAX=65536};

// NEW VARIABLE TYPES

typedef struct {
  unsigned int first;
  unsigned int isample;
  unsigned int Nrep;
  unsigned int NUM_WALKERS_PER_SITE;
  unsigned int Lv;
  unsigned int t;
  unsigned long long seed;
  double beta;
} icds;

typedef unsigned long long int MYWORD;
#define BITSINMYWORD 64

#ifdef MAIN
typedef struct{uint4 vec[NR];} s_time;
typedef struct{uint32_t my_key[NR];} s_keys;
#else
// CUDA types in C
typedef struct{uint32_t x,y,z,w;} uint4;

typedef struct{uint4 vec[NR];} s_time;
typedef struct{uint32_t my_key[NR];} s_keys;
#endif

#define NDATULL 2
#define NDATINT 6
#define NDATRAND (2+NR)

typedef struct
{
  unsigned long long int maxtime,
    now;
  int ntimes,
    nterm,
    nc,
    nr,
    l,
    flag;
  randint seed_J,
    seed_u,
    seed_MC[NR];
} s_data;

typedef struct{
  uint4 umbrales[128];
  uint4 prebusqueda[NPREBUSQUEDAS>>3];
} s_lut_heat_bath;

struct Vicini {
    MYWORD J[MAXNSAMPLES][DEGREE];
    int neig[DEGREE];
} ;

// VARIABLES
#ifdef MAIN
int list_samples[MAXNSAMPLES];
unsigned long long int *list_times;

int x_p[Lx],x_m[Lx],y_p[Ly],y_m[Ly],z_p[Lz],z_m[Lz];

Vicini viciniB[HALF_MSC_V], viciniN[HALF_MSC_V];
unsigned int bianchi_index[HALF_MSC_V], neri_index[HALF_MSC_V];
unsigned char bianchi_rotate[HALF_MSC_V], neri_rotate[HALF_MSC_V];
int side[DIM], off[DIM+1]; 

int numThreadsPerBlock, numBlocks;

s_data data;
int write_seeds;
double betas[NBETAS];
int countJ0[MAXNSAMPLES];

s_lut_heat_bath h_LUT[NBETAS];

s_keys seed_keys;
s_time s_time_and_entropy;
randint seeds_J[MAXNSAMPLES], seeds_C[MAXNSAMPLES*NR];
s_aleatorio_HQ_64bits random_u, random_J;
s_aleatorio_HQ_64bits random_PRC[NR], random_PRC_C[MAXNSAMPLES][NR];
s_xoshiro256pp random_xoshiro256pp[NR], random_xoshiro256pp_C[MAXNSAMPLES][NR];
s_aleatorio_HQ_64bits write_random_PRC_C[MAXNSAMPLES][NR];
s_xoshiro256pp write_random_xoshiro256pp_C[MAXNSAMPLES][NR];

char uu[MAXNSAMPLES][NR][NBETAS][V], Jx[MAXNSAMPLES][V], Jy[MAXNSAMPLES][V], Jz[MAXNSAMPLES][V];
MYWORD h_MSC_u_even[MAXNSAMPLES][HALF_MSC_VNRNBETAS], h_MSC_u_odd[MAXNSAMPLES][HALF_MSC_VNRNBETAS];
MYWORD **h_J[MAXNSAMPLES], **h_spin[MAXNSAMPLES];
int **h_neig[MAXNSAMPLES];
int *h_Ener;
unsigned long long ener_C[MAXNSAMPLES][NRNBETAS];

uint8_t which_clon_this_beta[MAXNSAMPLES][NR][NBETAS];
uint8_t which_beta_this_clon[MAXNSAMPLES][NR][NBETAS];

double average_Energy[MAXNSAMPLES][NBETAS], std_Energy[MAXNSAMPLES][NBETAS];
double average_Energy_C[MAXNSAMPLES][NBETAS], std_Energy_C[MAXNSAMPLES][NBETAS];

int sum[MAX];

MYWORD *CPUwalker0[MAXNSAMPLES], *CPUwalker1[MAXNSAMPLES];
MYWORD *write_CPUwalker0[MAXNSAMPLES], *write_CPUwalker1[MAXNSAMPLES];

uint64_t threshold[DIM];

icds *p2icds[MAXNSAMPLES];

// DEVICE VARIABLES
#define DEVICONST __device__ __constant__
#define DEVIBASE __device__

DEVICONST int d_deltaE, d_deltaN, d_deltaU, d_deltaO, d_deltaS, d_deltaD;
unsigned int *d_bianchi_index, *d_neri_index;
unsigned char *d_bianchi_rotate, *d_neri_rotate;

__managed__ MYWORD **d_J[MAXNSAMPLES];
__managed__ int **d_neig[MAXNSAMPLES];
__managed__ MYWORD  **d_spin[MAXNSAMPLES];

MYWORD **h_walker0[MAXNSAMPLES];
MYWORD **h_walker1[MAXNSAMPLES];
__managed__ MYWORD  **d_walker0[MAXNSAMPLES];
__managed__ MYWORD  **d_walker1[MAXNSAMPLES];

__managed__ unsigned int *d_whichclone[MAXNSAMPLES];

s_lut_heat_bath *dev_lut_heat_bath;
uint64_t *rand_wht_h, *rand_blk_h;

int *d_Ener;
int *d_sum;

#else
extern int list_samples[];
extern unsigned long long int *list_times;

extern int x_p[],x_m[],y_p[],y_m[],z_p[],z_m[];

extern Vicini viciniB[HALF_MSC_V], viciniN[HALF_MSC_V];
extern unsigned int bianchi_index[], neri_index[];
extern unsigned char bianchi_rotate[], neri_rotate[];
extern int side[], off[]; 

extern int numThreadsPerBlock, numBlocks;

extern s_data data;
extern int write_seeds;
extern double betas[];
extern int countJ0[];

extern s_lut_heat_bath h_LUT[];

extern s_keys seed_keys;
extern randint seeds_J[], seeds_C[];
extern s_aleatorio_HQ_64bits random_u, random_J;
extern s_aleatorio_HQ_64bits random_PRC[], random_PRC_C[][NR];
extern s_xoshiro256pp random_xoshiro256pp[NR], random_xoshiro256pp_C[][NR];
extern s_aleatorio_HQ_64bits write_random_PRC_C[][NR];
extern s_xoshiro256pp write_random_xoshiro256pp_C[][NR];


extern char uu[][NR][NBETAS][V], Jx[][V], Jy[][V], Jz[][V];
extern MYWORD h_MSC_u_even[][HALF_MSC_VNRNBETAS], h_MSC_u_odd[][HALF_MSC_VNRNBETAS];
extern int *h_Ener;
extern unsigned long long ener_C[MAXNSAMPLES][NRNBETAS];

extern uint8_t which_clon_this_beta[][NR][NBETAS];
extern uint8_t which_beta_this_clon[][NR][NBETAS];

extern double average_Energy[MAXNSAMPLES][NBETAS], std_Energy[MAXNSAMPLES][NBETAS];
extern double average_Energy_C[MAXNSAMPLES][NBETAS], std_Energy_C[MAXNSAMPLES][NBETAS];

extern int sum[];


extern MYWORD *CPUwalker0[], *CPUwalker1[];
extern MYWORD *write_CPUwalker0[], *write_CPUwalker1[];

extern uint64_t threshold[];

extern icds *p2icds[];

#endif

// FUNCTIONS PROTOTYPES

// IO:
void write_conf(unsigned long long,int);
void read_conf(int);
void Janus_packing_for_write(int, int);
void Janus_unpacking_for_read(int, int);
void write_walkers(unsigned long long, int, int);
void read_walkers(int);
void write_measures(unsigned long long, int);
void read_input(const char *);
void check_data(s_data);
void print_data(FILE *, s_data *);
void read_betas(const char *);
void read_lut(const char *, const char *);
void read_list_times(const char *);
void read_list_samples(const char *, int);
void create_sample_path(int, int);
void check_and_prepare_simulation_backup(int, int);
int get_seeds(int);
void write_seeds_in_file(int);
void backup(int);

void stop(unsigned long long int *);
void ignore_sigterm(void);
void handle_sigterm(void);
void stop_run(int);
void create_error(const char *format, ...);
void create_running(void);
void renew_running(void);
void delete_running(void);
void writelog(FILE *, time_t);
void print_help(const char *);
void print_and_exit(const char *format, ...);

#ifdef JANUS
int compact_dump(char*,FILE*,int);
void write_conf_Js(int);
#endif

// TIMING
void measure_time(int, unsigned long long);
int check_time(unsigned long long, unsigned long long, unsigned long long);

// INIT:
void Init(int);
void Init_Random(int);
void Init_neighbours(void);
void Set_MSC_neigh(unsigned int *, int,int, int);
int punto(int *, int *);
void coordinate(int, int *,int *);
void sp(int *, int * , int *, int);
void Init_MSC_neighbours(void);
void Init_MSC_index_and_rotate(void);
void Init_Binario(void);

void Init_u(int);
void packing_u(int);
void unpacking_u(int);
void check_packing(int);

void Init_J(int);
void packing_J(int);

void calculate_blocks(void);

void init_tempering(int);

void Init_threads_for_walkers(int);
void get_thresholds(uint64_t *, double);
void * InitWalkers(void *);

// MED:
void Measures(int);
void Measure_Energy_Clerics(int);
int calculate_scalar_energy(int, int, int);

// MAIN:
int Monte_Carlo(uint32_t, unsigned long long, int, int);
void Measure_Energy_Spins(int);

// DEVICE MEMORY
void refresh_walkers(int);
void CopyConstOnDevice(void);
void CopyDataOnDevice(int);
void inline copy_from_GPU_for_write(int);
void send_PT_state_to_GPU(int);
void FreeDevice(int);
void FreeHost(int);

// KERNELS
#ifdef MAIN
__global__ void d_gen_rndbits_cuda(int, int, int, s_time, s_keys, s_lut_heat_bath *, uint64_t *);

__global__ void d_computeEne(MYWORD *, MYWORD *, int *, int *, int, MYWORD **, int **, int,
			     int, int, unsigned int *, unsigned int *);

__global__ void d_oneMCstepBN_multisample(MYWORD **[],
					  const int,
					  MYWORD **[],
					  int **[],
					  int,
					  int,
					  int,
					  unsigned int **,
					  const unsigned char *,
					  uint64_t *,
					  MYWORD **[], MYWORD **[], int);


__global__ void d_Walk(const int , int **[],
			       int ,
			       int ,
			       int ,
			       int ,
			       const unsigned char *,
			       MYWORD **[],
			       MYWORD **[]);

#ifdef MASSIMO
__device__ void Sum(MYWORD * sum_bit3, MYWORD *sum_bit2, MYWORD * sum_bit1,
		    MYWORD num1_bit3, MYWORD num1_bit2, MYWORD num1_bit1,
		    MYWORD num2_bit3, MYWORD num2_bit2, MYWORD num2_bit1);

#endif
#endif
void **mmcuda(void ***, int , int , int , int);

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <stdarg.h>
#define NBR 64
#include "random.h"


// Como no nos fiamos del todo de _philox_4_x32_10, intentamos mejorar la calidad añadiendo
// en cada llamada del kernel un cierto número de bits de entropia, generados en la CPU mediante
// un generador de Parisi-Rapuano+Congruencial (ambos de 64 bits) sumados a un xoshiro256++.
//
// La cuenta de bits es como sigue.
//
// El numero de bits del contador de _philox_4x32_10 es 128, divididos en 4 palabras de 32.
//
// Con el método que queremos utilizar hacen falta 3 llamadas a _philox_4x32_10 dentro del kernel.
// Es decir, hay que reservar 2 bits para el contador interno, que se actualiza en el kernel.
//
// Reservamos 32 bits para la iteración de Metropolis, mas el bit de paridad hacen 33. (Es un factor 16
// de seguridad respecto a lo que Isi cree que necesita).
//
// La hebra (tid) necesita 21 bits: 86 betas * 2^{14} palabras de 32 bits=2^{20.42...}
//
// El jobID precisa 14 bits distribuidos como: 6 replicas (3 bits) + 1280 samples (11 bits)
//
// Por tanto, nos quedan 128-(2+32+1+14+21)=58 bits sin utilizar en el contador que podemos emplear
// para inyectar entropía desde los generadores de la CPU en la secuencia aleatoria de la GPU. Vamos a distribuir
// esos 47 bits en cada una de las 4 palabras del contador.


//#define NUMBETAS 86
#ifndef NUMBETAS
#error "NUMBETAS no definido"
#endif

#define NPREBUSQUEDAS (1<<NUMBITSPREBUSQUEDAS)//Tiene que ser una potencia de dos
#define NUMBITS_THREADS_PER_BLOCK 10
#define THREADSPERBLOCK (1<<NUMBITS_THREADS_PER_BLOCK)
#define NUMBLOCKS (16*NUMBETAS) 

// Flexibilidad cero con estos define
#define Lx 16
#define Ly 16
#define Lz 4096
#define MSC_Lz (Lz>>6)
#define MSC_V (Lx*Ly*MSC_Lz)
#define V (Lx*Ly*Lz)

//Esta es la información para nuestro "heat-bath" de bits aleatorios 
typedef struct{
  uint4 umbrales[128];  //256 umbrales con precision de 64 bits
  uint4 prebusqueda[NPREBUSQUEDAS>>3]; //Para disminuir el número de iteracciones en la bisección.
} s_lut_heat_bath;

#if(THREADSPERBLOCK<(128+(NPREBUSQUEDAS>>3)))
#error "El esquema de paso a _shared_ no funciona, aumenta THREADSPERBLOCK"
#endif



//Rutinas de inicializacion:
void Init(void);
void acera_acumuladores(void);
//Rutinas de comprobación:
void comprueba_resultados(long long int, uint64_t, uint32_t);
//Rutinas de IO:
void lee_betas(char *);
void lee_lut(char *);
void print_and_exit(const char *format, ...);
void obten_crc(uint32_t *);


//Rutinas para el "run_time" de CUDA (caja negra tomada del "cuda_by_example")

static void HandleError( cudaError_t err,
                         const char *file,
                         int line ) {
    if (err != cudaSuccess) {
        printf( "%s in %s at line %d\n", cudaGetErrorString( err ),
                file, line );
        exit( EXIT_FAILURE );
    }
}
#define HANDLE_ERROR( err ) (HandleError( err, __FILE__, __LINE__ ))


#define HANDLE_NULL( a ) {if (a == NULL) { \
                            printf( "Host memory failed in %s at line %d\n", \
                                    __FILE__, __LINE__ ); \
                            exit( EXIT_FAILURE );}}


//Variables globales

#ifdef MAIN

s_lut_heat_bath lut_heat_bath[NUMBETAS];
double betas[NUMBETAS];
int cual_beta_este_clon[NUMBETAS];
int cual_clon_esta_beta[NUMBETAS];
uint4 suma_global[NUMBETAS*MSC_V];

//Punteros de la GPU

uint4 * dev_suma_global; //Longitud: NUMBETAS*MSC_V
s_lut_heat_bath * dev_lut_heat_bath; //Longitud: NUMBETAS
__device__ __constant__ int32_t dev_cual_beta_este_clon[NUMBETAS];

#else

extern s_lut_heat_bath lut_heat_bath[NUMBETAS];
extern double betas[NUMBETAS];
extern int cual_beta_este_clon[NUMBETAS];
extern int cual_clon_esta_beta[NUMBETAS];
extern uint4 suma_global[NUMBETAS*MSC_V];

#endif

#define MAIN
#include "bits-random.h"

void  alloca_memoria_GPU(void)
{

  size_t total;
  
  total=NUMBETAS*MSC_V*sizeof(uint4);
  HANDLE_ERROR( cudaMalloc((void **)&dev_suma_global,total));

  total=NUMBETAS*sizeof(s_lut_heat_bath);
  HANDLE_ERROR( cudaMalloc((void **)&dev_lut_heat_bath,total));

}


void  transfiere_memoria_CPU_a_GPU(void)
{

  size_t total;

  total=NUMBETAS*MSC_V*sizeof(uint4);
  HANDLE_ERROR( cudaMemcpy( dev_suma_global, suma_global, total , cudaMemcpyHostToDevice) );

  total=NUMBETAS*sizeof(s_lut_heat_bath);
  HANDLE_ERROR( cudaMemcpy( dev_lut_heat_bath, lut_heat_bath, total , cudaMemcpyHostToDevice) );

}

void libera_memoria_GPU(void)
{

  HANDLE_ERROR( cudaFree( dev_suma_global ));
  HANDLE_ERROR( cudaFree( dev_lut_heat_bath ));

}


//Parallel Tempering, in the context of this program, is only good for debug.

void Parallel_Tempering(s_xoshiro256pp * p)
{
  int ibeta,iclon,temp,cambiar;
  //double pt,exppt; //Variables needed by true Parallel Tempering


  for(ibeta=0;ibeta<(NUMBETAS-1);ibeta++){

    /*
    pt=-(betas[ibeta+1]-betas[ibeta])
      (ener_para_el_PT[cual_clon_esta_beta[ibeta+1]]-ener_para_el_PT[cual_clon_esta_beta[ibeta]]);
    */

    _actualiza_xoshiro256pp(p[0]);
    
    cambiar=(int)(p[0].final>>63);

    /* True Parallel Tempering*/
    /*
    cambiar=0;
    if (pt>=0){
      aceptanciaf[ibeta]++;
      cambiar=1;
    }else{
      exppt=exp(pt);      
      aceptanciaf[ibeta]+=exppt;
      if((FNORM*p[0].final) < exppt )
	cambiar=1;
    }
    */

    if(cambiar){ // aceptamos el cambio
      temp=cual_clon_esta_beta[ibeta];
      cual_clon_esta_beta[ibeta]=cual_clon_esta_beta[ibeta+1];
      cual_clon_esta_beta[ibeta+1]=temp;      
      //aceptancia[ibeta]++; //Descomentar en el verdadero programa
    }
  }

  //Actualizamos la permutacion inversa
  for (ibeta=0;ibeta<NUMBETAS;ibeta++){
    iclon=cual_clon_esta_beta[ibeta];    
    cual_beta_este_clon[iclon]=ibeta; 
  }

}


//ATENCION: en CUDA debemos usar mulhi en vez del multiplicador de 64 bits
//Defines for philox_4x32_10 implementation (homemade)
#ifndef PHILOX_M4x32_0
#define PHILOX_M4x32_0 ((uint32_t)0xD2511F53)
#endif
#ifndef PHILOX_M4x32_1
#define PHILOX_M4x32_1 ((uint32_t)0xCD9E8D57)
#endif

#ifndef PHILOX_W32_0
#define PHILOX_W32_0 ((uint32_t)0x9E3779B9)
#endif
#ifndef PHILOX_W32_1
#define PHILOX_W32_1 ((uint32_t)0xBB67AE85)
#endif


#define _actualiza_clave {\
    key[0]+=PHILOX_W32_0; \
    key[1]+=PHILOX_W32_1;}

#define _actualiza_estado {\
  lo0=PHILOX_M4x32_0*v[0];			\
  hi0=__umulhi(PHILOX_M4x32_0,v[0]);		\
  lo1=PHILOX_M4x32_1*v[2];			\
  hi1=__umulhi(PHILOX_M4x32_1,v[2]);		\
  v[0]=hi1^v[1]^key[0];				\
  v[1]=lo1;					\
  v[2]=hi0^v[3]^key[1];				\
  v[3]=lo0;}

#define _philox_4x32_10 {\
  _actualiza_estado; \
  _actualiza_clave;   \
  _actualiza_estado; \
  _actualiza_clave;   \
  _actualiza_estado; \
  _actualiza_clave;   \
  _actualiza_estado; \
  _actualiza_clave;   \
  _actualiza_estado; \
  _actualiza_clave;   \
  _actualiza_estado; \
  _actualiza_clave;   \
  _actualiza_estado; \
  _actualiza_clave;   \
  _actualiza_estado; \
  _actualiza_clave;   \
  _actualiza_estado; \
  _actualiza_clave;   \
  _actualiza_estado;}

#define _obten_aleatorio_contador {					\
    uint32_t hi0,hi1,lo0,lo1;						\
    uint32_t key[2];							\
    key[0]=tid; key[1]=useed;						\
    v[0]=tiempo_y_entropia.x;						\
    v[1]=tiempo_y_entropia.y;						\
    v[2]=tiempo_y_entropia.z;						\
    v[3]=tiempo_y_entropia.w;						\
    _philox_4x32_10;							\
    tiempo_y_entropia.x++; /* contador interno*/			\
    /* lsb: iteracción de xoroshiro128++ deconstruida por VMM */	\
    hi1=v[0]^v[2];							\
    lo1=v[1]^v[3];							\
    lsb[0]=((v[1]<<17)|(v[0]>>15))^hi1^((hi1<<21)|lo1>>11);		\
    lsb[1]=((v[0]<<17)|(v[1]>>15))^lo1^(lo1<<21);			\
    lsb[2]=(hi1<<28)|(lo1>>4);						\
    lsb[3]=(lo1<<28)|(hi1>>4);}


#define _biseccion(Rmsb,Rlsb) {						\
    const uint32_t msb=Rmsb>>(32-NUMBITSPREBUSQUEDAS);			\
    uint32_t min,max,decision;						\
    decision=(uint32_t) (~(((Rmsb==dev_umbrales[254].x)?(Rlsb<dev_umbrales[254].y):(Rmsb<dev_umbrales[254].x))-1)); \
    min=dev_prebusqueda[msb]&255;					\
    max=dev_prebusqueda[msb]>>8;					\
    min=(decision&min)|((~decision)&255);				\
    max=(decision&max)|((~decision)&255);				\
    while((max-min)>1){							\
      n=(max+min)>>1;							\
      decision=(uint32_t) (~(((Rmsb==dev_umbrales[n].x)?(Rlsb<dev_umbrales[n].y):(Rmsb<dev_umbrales[n].x))-1)); \
      max=(decision&n)|((~decision)&max);				\
      min=(decision&min)|((~decision)&n);				\
    }									\
  decision=(uint32_t) (~(((Rmsb==dev_umbrales[min].x)?(Rlsb<dev_umbrales[min].y):(Rmsb<dev_umbrales[min].x))-1)); \
  n=(decision&min)|(~decision&max);}

#define _obten_palabra(bb) {				\
    uint32_t v[4],lsb[4];				\
    _obten_aleatorio_contador;				\
    const uint32_t rot=(lsb[0]+v[1]+lsb[2]+v[3])>>27;	\
    _biseccion(v[0],lsb[0]);				\
    bb=n;						\
    _biseccion(v[1],lsb[1]);				\
    bb|=n<<8;						\
    _biseccion(v[2],lsb[2]);				\
    bb|=n<<16;						\
    _biseccion(v[3],lsb[3]);				\
    bb|=n<<24;						\
    bb=(bb<<rot)|(bb>>(32-rot));}


__global__ void obten_bits(uint4 tiempo_y_entropia,uint32_t useed,
			   s_lut_heat_bath * dev_lut_heat_bath,
			   uint4 * dev_suma_global)
{

  int tid=threadIdx.x|(blockIdx.x<<NUMBITS_THREADS_PER_BLOCK);
  uint4 suma_local;
  int ic,ibeta,site,j;
  uint32_t n,b4,b8,b12;
  __shared__ uint2 dev_umbrales[256];
  __shared__ uint32_t dev_prebusqueda[NPREBUSQUEDAS];

  

  //Esto depende crucialmente de que MSC_V=1<<14;
  ic=tid>>14;
  ibeta=dev_cual_beta_este_clon[ic]; //Lectura de constant memory

  //Ya sabemos quien es nuestra beta. Tenemos que copiar las LUT a la memoria shared.

  if(threadIdx.x<(128+(NPREBUSQUEDAS>>3))){
    uint4 palabra;      
    if(threadIdx.x<128){//Estas hebras reunen umbrales
      palabra=dev_lut_heat_bath[ibeta].umbrales[threadIdx.x];
      j=threadIdx.x<<1;
      dev_umbrales[j].x=palabra.x;
      dev_umbrales[j].y=palabra.y;
      j++;
      dev_umbrales[j].x=palabra.z;
      dev_umbrales[j].y=palabra.w;      
    }else{
      j=threadIdx.x-128;	
      palabra=dev_lut_heat_bath[ibeta].prebusqueda[j];
      j<<=3;
      dev_prebusqueda[j++]=palabra.x&0xffffU;
      dev_prebusqueda[j++]=palabra.x>>16;
      dev_prebusqueda[j++]=palabra.y&0xffffU;
      dev_prebusqueda[j++]=palabra.y>>16;
      dev_prebusqueda[j++]=palabra.z&0xffffU;
      dev_prebusqueda[j++]=palabra.z>>16;
      dev_prebusqueda[j++]=palabra.w&0xffffU;
      dev_prebusqueda[j]=palabra.w>>16;
    }
  }

  __syncthreads(); //Las threads deben sincronizarse todas aquí
  
  //Empieza el cálculo
  
  site=tid&(MSC_V-1);//MSC_V debe ser potencia de dos!
  tiempo_y_entropia.z+=tid;
  suma_local=dev_suma_global[(ibeta<<14)|site];

  _obten_palabra(b4);
  _obten_palabra(b8);
  _obten_palabra(b12);
  b8=b4&b8;
  b12=b8&b12;
  
  suma_local.x+=__popc(b4);
  suma_local.y+=__popc(b8);
  suma_local.z+=__popc(b12);

  dev_suma_global[(ibeta<<14)|site]=suma_local;
  
}
#undef _obten_palabra
#undef _obten_aleatorio_contador
#undef _philox_4x32_10
#undef _actualiza_clave
#undef _actualiza_estado
#undef PHILOX_M4x32_0
#undef PHILOX_M4x32_1
#undef PHILOX_W32_0
#undef PHILOX_W32_1



//MAIN: (por comodidad, la colocación de los bits del contador la hacemos en un macro aparte)
//Los 32 bits de tiempo se distribuyen en (15+bit de paridad) en .x y 17 en .y
//En .x tenemos 2 bits de contador interno y (15+1) bits de contador externo. Son 18 bits.
//Queda sitio para 14 bits de entropia.
//En .y tenemos 17 bits de contador externo, nos queda sitio para 15 bits de entropia.
//En .z tenemos reservados 21 bits para el índice de thread, tid. Caben 11 bits de entropía.
//En .w tenemos reservados 14 bits para el índice de jobID. Nos queda sitio para 18 bits de entropía
#define _coloca_bits {                                                  \
    _actualiza_xoshiro256pp(aleatorio_xoshiro256pp);                    \
    _actualiza_aleatorio_HQ_escalar(aleatorio_PRC);                     \
    entropia=aleatorio_PRC.final+aleatorio_xoshiro256pp.final;          \
    entropia>>=6; /*Dejamos los 58 bits */                              \
    temporal=iter>>15;                                                  \
    tiempo_y_entropia.y=(uint32_t) temporal;                            \
    temporal=(iter<<49)>>46; /* Colocamos 15 bits a partir del bit 3*/  \
    temporal|=(paridad&1ULL)<<2;                                        \
    tiempo_y_entropia.x=(uint32_t) temporal;                            \
    temporal=(entropia&16383ULL)<<18;					\
    tiempo_y_entropia.x|=(uint32_t) temporal;                           \
    temporal=((entropia>>14)&32767ULL)<<17;                             \
    tiempo_y_entropia.y|=(uint32_t) temporal;                           \
    temporal=((entropia>>29)&2047ULL)<<21;				\
    tiempo_y_entropia.z=(uint32_t) temporal;                            \
    temporal=(entropia>>40)<<14;                                        \
    tiempo_y_entropia.w=((uint32_t) temporal)|mi_jobID;}


int main(int argc, char **argv)
{
  unsigned long long semilla_leida;
  uint64_t semilla,iter,entropia,paridad;
  uint32_t mi_jobID,semilla_key;
  uint4 tiempo_y_entropia;
  uint64_t temporal;  
  s_aleatorio_HQ_64bits aleatorio_PRC;
  s_xoshiro256pp aleatorio_xoshiro256pp;
  long long int num_pasos,num_PT,PTfrec;
  int iMetropolis;
  int device,count_devices;
  char name_betas[1024];
  size_t total;
  cudaEvent_t start, stop;
  float elapsedTime;

  
  if(argc!=6)
    print_and_exit("Usage: %s semilla numPT Ptfrec name_betas device\n",argv[0]);

  sscanf(argv[1],"%llu",&semilla_leida);
  sscanf(argv[2],"%lld",&num_PT);
  sscanf(argv[3],"%lld",&PTfrec);
  num_pasos=num_PT*PTfrec;
  sscanf(argv[4],"%s",name_betas);
  sscanf(argv[5],"%d",&device);
  HANDLE_ERROR( cudaGetDeviceCount( &count_devices) );
  if((device<0)||(device>=count_devices))
    print_and_exit("device=%d count_devices=%d\n",device,count_devices);
  
  HANDLE_ERROR( cudaSetDevice(device) );
  
  lee_betas(name_betas);
  lee_lut(name_betas);


  semilla=comprueba_semilla((uint64_t) semilla_leida);
  printf("Comprobando el generador de la CPU con semilla %llu\n", (unsigned long long) semilla);
  
  Inicia_generadores_CPU(&aleatorio_PRC,&aleatorio_xoshiro256pp,&semilla_key,semilla);

  printf("semilla_key=%u (binario:%u%u%u%u...)\n",semilla_key,(semilla_key>>31)&1U,
	 (semilla_key>>30)&1U,(semilla_key>>29)&1U,(semilla_key>>28)&1U);

  Init();

  //Alloca memoria para el dispositivo

  alloca_memoria_GPU();
  transfiere_memoria_CPU_a_GPU();

  
  //Generamos los bits random en la GPU

  
  mi_jobID=(12U<<10)+709U; //Por ejemplo. Se trata de un número representable con 14 bits.
  iter=0;
  total=NUMBETAS*sizeof(int32_t);

  //Comenzamos la cuenta de tiempos
  HANDLE_ERROR( cudaEventCreate( &start ) );
  HANDLE_ERROR( cudaEventCreate( &stop ) );
  HANDLE_ERROR( cudaEventRecord( start, 0 ) );
  
  for(iter=0;iter<num_pasos;){


    //Copia de la permutación del PT a la GPU
    HANDLE_ERROR( cudaMemcpyToSymbol ( dev_cual_beta_este_clon, cual_beta_este_clon, total, (size_t) 0, cudaMemcpyHostToDevice ) );

    
    //Así se mandarán las cosas
    for(iMetropolis=0;iMetropolis<PTfrec;iMetropolis++){
      paridad=0;
      _coloca_bits;   
      //comprueba(tiempo_y_entropia,iter,paridad,mi_jobID,entropia);
      obten_bits<<<NUMBLOCKS,THREADSPERBLOCK>>>(tiempo_y_entropia,
						semilla_key,
						dev_lut_heat_bath,
						dev_suma_global);
    
      paridad=1;   
      _coloca_bits;
      //comprueba(tiempo_y_entropia,iter,paridad,mi_jobID,entropia);
      obten_bits<<<NUMBLOCKS,THREADSPERBLOCK>>>(tiempo_y_entropia,
						semilla_key,
						dev_lut_heat_bath,
						dev_suma_global);

      iter++;
    }
    
    Parallel_Tempering(&aleatorio_xoshiro256pp);
  }


  HANDLE_ERROR( cudaEventRecord( stop, 0 ) );
  HANDLE_ERROR( cudaEventSynchronize( stop ) );
  HANDLE_ERROR( cudaEventElapsedTime( &elapsedTime, start, stop ) );
  printf("Dispositivo numero: %d\n",device);
  printf( "Tiempo empleado por la GPU en obtener bits random: %3.1f ms\n", elapsedTime );
  double perspinupdate=((double) elapsedTime)/((double) V*NUMBETAS*num_pasos);
  perspinupdate*=1e9;
  printf("Tiempo consumido en bits random por update: %lf ps\n",perspinupdate);
 
  HANDLE_ERROR( cudaEventDestroy( start ) );
  HANDLE_ERROR( cudaEventDestroy( stop ) );

  
  //Leemos los resultados de la GPU
  total=MSC_V*NUMBETAS*sizeof(uint4);
  HANDLE_ERROR( cudaMemcpy(suma_global, dev_suma_global, total, cudaMemcpyDeviceToHost));

  comprueba_resultados(num_pasos,semilla,semilla_key);

  //Liberamos memoria
  libera_memoria_GPU();

  return 0;
  
}

#undef _coloca_bits

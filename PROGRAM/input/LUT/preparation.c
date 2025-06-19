#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <quadmath.h>
#include <stdarg.h>
#include <stdint.h>

//#define NUMBETAS 86
#ifndef NUMBETAS
#error "NUMBETAS no definido"
#endif


#define NPREBUSQUEDAS (1<<NUMBITSPREBUSQUEDAS)//Tiene que ser una potencia de dos
#define NUMBITS_THREADS_PER_BLOCK 10
#define THREADSPERBLOCK (1<<NUMBITS_THREADS_PER_BLOCK)


typedef struct{uint32_t x,y,z,w;} uint4;
typedef struct{uint32_t x,y;} uint2;

typedef struct{
  uint4 umbrales[128];
  uint4 prebusqueda[NPREBUSQUEDAS>>3];
} s_lut_heat_bath;


#if(THREADSPERBLOCK<(128+(NPREBUSQUEDAS>>3)))
#error "El esquema de paso a _shared_ no funciona, aumenta THREADSPERBLOCK"
#endif


s_lut_heat_bath lut_heat_bath[NUMBETAS];
double betas[NUMBETAS];


void print_and_exit(const char *format, ...);
void prepara_LUT(double, s_lut_heat_bath *);
void lee_betas(char *);
int desempaqueta_y_busca_por_biseccion(int);
void escribe_luts(void);


int main(int argc,char **argv)
{
  int ibeta,i,j,min,max,dmax,d,maxiter;
  char name_betas[1024];
  union{
    uint4 palabra;
    uint32_t vec[4];
  }paquete;

  
  if(argc!=2)
    print_and_exit("Usage: %s name_betas\n",argv[0]);
  sscanf(argv[1],"%s",name_betas);
  lee_betas(name_betas);

  for(ibeta=0;ibeta<NUMBETAS;ibeta++)
    prepara_LUT(betas[ibeta],lut_heat_bath+ibeta);

  printf("\n"); printf("Eficacia de prebusqueda:\n");
  for(ibeta=0;ibeta<NUMBETAS;ibeta++){
    dmax=-1;

    for(i=0;i<NPREBUSQUEDAS;i++){
      paquete.palabra=lut_heat_bath[ibeta].prebusqueda[i>>3];
      j=(i&1)<<4;
      min=(paquete.vec[(i&7)>>1]>>j)&255;
      max=(paquete.vec[(i&7)>>1]>>(j+8))&255;
      d=max-min;
      if(d<0){
	print_and_exit("Fallo de ordenacion beta=%.14g msb=%d\n",
		       betas[ibeta],i);
      }
      dmax=(dmax<d)?d:dmax;      
    }
    maxiter=desempaqueta_y_busca_por_biseccion(ibeta);
    printf("beta=%.14g dmax=%d, maxiter=%d\n",betas[ibeta],d,maxiter);
  }

  escribe_luts();
  return 0;
}


int desempaqueta_y_busca_por_biseccion(int ibeta)
{
  
  int j,msb;
  uint32_t Rmsb,Rlsb,min,max,n,decision;
  //Estas dos serán _shared_ en la GPU
  uint2 dev_umbrales[256];
  uint32_t dev_prebusqueda[NPREBUSQUEDAS];
  //Estas variables no hacen falta en la GPU
  int maxiter,niter;
  uint32_t n_debe_ser;
  uint64_t Rmin,Rmax;
  struct {int x;} threadIDx;
  
  //Primero desempaquetar, emulemos lo que ocurre en la GPU;
  for(threadIDx.x=0;threadIDx.x<THREADSPERBLOCK;threadIDx.x++){

    if(threadIDx.x<(128+(NPREBUSQUEDAS>>3))){
      uint4 palabra;      
      if(threadIDx.x<128){//Estas hebras reunen umbrales
	palabra=lut_heat_bath[ibeta].umbrales[threadIDx.x];
	j=threadIDx.x<<1;
	dev_umbrales[j].x=palabra.x;
	dev_umbrales[j].y=palabra.y;
	j++;
	dev_umbrales[j].x=palabra.z;
	dev_umbrales[j].y=palabra.w;      
      }else{
	j=threadIDx.x-128;	
	palabra=lut_heat_bath[ibeta].prebusqueda[j];
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
  }
  // __syncthreads(); Las threads deben sincronizarse todas aquí

  //Hacemos una búsqueda exhaustiva en cada intervalo de probabilidad.
  maxiter=-1;
  
  for(n_debe_ser=0;n_debe_ser<256;n_debe_ser++){
    if(n_debe_ser==0){
      Rmin=0;
      Rmax=(((uint64_t) dev_umbrales[n_debe_ser].x)<<32)|((uint64_t) dev_umbrales[n_debe_ser].y);
      Rmax--;
    }else{
      if(n_debe_ser==255){
	Rmax=~0ULL;
	Rmin=(((uint64_t) dev_umbrales[n_debe_ser-1].x)<<32)|((uint64_t) dev_umbrales[n_debe_ser-1].y);
      }else{
	Rmin=(((uint64_t) dev_umbrales[n_debe_ser-1].x)<<32)|((uint64_t) dev_umbrales[n_debe_ser-1].y);
	Rmax=(((uint64_t) dev_umbrales[n_debe_ser].x)<<32)|((uint64_t) dev_umbrales[n_debe_ser].y);
	Rmax--;
      }
    }
    
    //Empezamos con Rmin:
    Rmsb=(uint32_t)(Rmin>>32); Rlsb=(uint32_t)(Rmin&0xffffffffU);

    //Esta es la bisección
    msb=Rmsb>>(32-NUMBITSPREBUSQUEDAS);
    decision=(uint32_t) (~(((Rmsb==dev_umbrales[254].x)?(Rlsb<dev_umbrales[254].y):
			    (Rmsb<dev_umbrales[254].x))-1));

    min=dev_prebusqueda[msb]&255;
    max=dev_prebusqueda[msb]>>8;
    min=(decision&min)|((~decision)&255);
    max=(decision&max)|((~decision)&255);
    niter=0;
    while((max-min)>1){
      n=(max+min)>>1;
      decision=(uint32_t) (~(((Rmsb==dev_umbrales[n].x)?(Rlsb<dev_umbrales[n].y):
			      (Rmsb<dev_umbrales[n].x))-1));
      max=(decision&n)|((~decision)&max);
      min=(decision&min)|((~decision)&n);
      niter++;
    }
    niter++;
    decision=(uint32_t) (~(((Rmsb==dev_umbrales[min].x)?(Rlsb<dev_umbrales[min].y):
			    (Rmsb<dev_umbrales[min].x))-1));
    n=(decision&min)|(~decision&max);
    if(n!=n_debe_ser){
      printf("beta=%.14g, la pifiamos en la bisección\n",betas[ibeta]);
      printf("Rmin=%llu Rmsb=%u Rlsb=%u\n",(unsigned long long) Rmin,Rmsb,Rlsb);
      print_and_exit("n_debe_ser=%u n=%u\n",n_debe_ser,n);
    }

    maxiter=(niter>maxiter)?niter:maxiter;

    //Seguimos con Rmax:
    Rmsb=(uint32_t)(Rmax>>32); Rlsb=(uint32_t)(Rmax&0xffffffffU);
    //Esta es la bisección
    msb=Rmsb>>(32-NUMBITSPREBUSQUEDAS);
        decision=(uint32_t) (~(((Rmsb==dev_umbrales[254].x)?(Rlsb<dev_umbrales[254].y):
			    (Rmsb<dev_umbrales[254].x))-1));

    min=dev_prebusqueda[msb]&255;
    max=dev_prebusqueda[msb]>>8;
    min=(decision&min)|((~decision)&255);
    max=(decision&max)|((~decision)&255);

    niter=0;
    while((max-min)>1){
      n=(max+min)>>1;
      decision=(uint32_t) (~(((Rmsb==dev_umbrales[n].x)?(Rlsb<dev_umbrales[n].y):
			      (Rmsb<dev_umbrales[n].x))-1));
      max=(decision&n)|((~decision)&max);
      min=(decision&min)|((~decision)&n);
      niter++;
    }
    decision=(uint32_t) (~(((Rmsb==dev_umbrales[min].x)?(Rlsb<dev_umbrales[min].y):
			    (Rmsb<dev_umbrales[min].x))-1));
    n=(decision&min)|(~decision&max);
    if(n!=n_debe_ser){
      printf("beta=%.14g, la pifiamos en la bisección\n",betas[ibeta]);
      printf("Rmax=%llu Rmsb=%u Rlsb=%u\n",(unsigned long long) Rmax,Rmsb,Rlsb);
      print_and_exit("n_debe_ser=%u n=%u\n",n_debe_ser,n);
    }

    maxiter=(niter>maxiter)?niter:maxiter;
  }
      
  return maxiter;
}

void lee_betas(char * name_betas)
{
  int ibeta;

  double dummy;
  FILE *fi_betas;
  

  
  if (NULL==(fi_betas=fopen(name_betas,"rt")))
    print_and_exit(" No existe el fichero %s.\n",name_betas);

  for (ibeta=0;ibeta<NUMBETAS;ibeta++)
    if (EOF==fscanf(fi_betas,"%lf ",&betas[ibeta]))
      print_and_exit("El fichero %s solo tiene %d lineas (necesita %d)\n",
		     name_betas,ibeta,NUMBETAS);
  // hay alguna beta de sobra en el fichero?
  if (EOF!=fscanf(fi_betas,"%lf",&dummy))
    fprintf(stderr,"Hay lineas de mas en el fichero %s\n",name_betas); 
  
  for(ibeta=0;ibeta<NUMBETAS-1;ibeta++){
    if (betas[ibeta]<betas[ibeta+1])
      print_and_exit("Error en fichero %s: betas no ordenadas de mayor a menor!\n",name_betas);
  }
  
}


#define TWOBRQ 18446744073709551616.Q //2^64

//Sea numbits[n] el popcount de n, con 0<= n < 256
//
// La probabilidad de encontrar "n" es:
//
// prob[n] = exp(-4*beta*numbits*numbits[n])*
//              (1-exp(-4*beta*numbits)^{8-numbits[n]);
//
// umbral[n]=\sum_{i=0}^n 2^{64}* prob[n].
//
// Sea R un random uniforme de 64 bits.
//
// Si umbral[n-1] <= R< umbral[n], elegiremos "n".
//
// Tanto los umbrales como R se trocean en sus 32 msb y lsb
// Para la prebusqueda utilizaremos los 8 bits más significativos de R. 

void prepara_LUT(double beta, s_lut_heat_bath * mi_lut)
{
  int num_bits,i,j;
  int histograma[9];
  uint32_t tabla_min[NPREBUSQUEDAS],tabla_max[NPREBUSQUEDAS];
  uint64_t msb,cola;
  uint64_t caso[9];
  uint64_t cumulativa[256];
  uint64_t ultimo_intervalo,temporal;
  __float128 qtemp,qbeta,qprob0,qprob1,probabilidad,total;
  __float128 qcaso[9];

  union{
    uint4 palabra;
    uint32_t vec[4];
  }paquete;

  qbeta=(__float128) beta;
  qprob0=(1.0q-expq(-4.0q*qbeta));
  qprob1=expq(-4.0q*qbeta);
  total=0;
  for(i=0;i<9;i++)
    histograma[i]=0;
  for(i=0;i<256;i++)
    histograma[__builtin_popcount(i)]++;
  
  for(i=0;i<9;i++){
    probabilidad=powq(qprob0,8-i)*powq(qprob1,i);
    qcaso[i]=probabilidad;
    total+=qcaso[i]*((__float128) histograma[i]);
    //printf("i=%d histo[%d]=%d total=%.14g\n",i,i,histograma[i],(double) total);
    qtemp=TWOBRQ*probabilidad;
    caso[i]=(uint64_t) qtemp;
  }

  
  cumulativa[0]=caso[0];
  total=qcaso[0];
  for(i=1;i<255;i++){
    num_bits=__builtin_popcount(i);
    total+=qcaso[num_bits];
    cumulativa[i]=cumulativa[i-1]+caso[num_bits];
  }
  cumulativa[255]=~0ULL;
  total+=qcaso[8];
  total-=1.0q;
  printf("beta=%.14g, defecto de normalizacion 128 bits: %.14g. ",
	 beta,(double) total);
  
  ultimo_intervalo=cumulativa[255]-cumulativa[254]+1;
  printf("Ultimo intervalo debe: %llu pero es %llu (cociente-1: %.14g)\n",
	 (unsigned long long) caso[8],
	 (unsigned long long) ultimo_intervalo,
	 -1.+((double) ultimo_intervalo)/((double) caso[8]));

  for(i=0;i<256;i+=2){    
    paquete.vec[0]=(uint32_t) (cumulativa[i]>>32);
    paquete.vec[1]=(uint32_t) (cumulativa[i]&0xffffffffULL);
    paquete.vec[2]=(uint32_t) (cumulativa[i+1]>>32);
    paquete.vec[3]=(uint32_t) (cumulativa[i+1]&0xffffffffULL);
    mi_lut[0].umbrales[i>>1]=paquete.palabra;    
  }
  
  //Comprobamos el empaquetamiento:
  for(i=0;i<256;i++){
    paquete.palabra=mi_lut[0].umbrales[i>>1];    
    temporal=(((uint64_t) paquete.vec[(i&1)<<1])<<32)|
      ((uint64_t) paquete.vec[((i&1)<<1)|1]);
    if(temporal!=cumulativa[i])
      print_and_exit("Fallo de empaquetamiento, beta=%.14 g, i=%d\n",
		     beta,i);
  }
    
  //Ahora hacemos la prebusqueda

  cola=0ULL;
  for(i=0;i<(64-NUMBITSPREBUSQUEDAS);i++)
    cola|=(1ULL<<i);

  if(cola!=((1LL<<(64-NUMBITSPREBUSQUEDAS))-1LL))
    print_and_exit("cola=%llu (debe ser %ll)\n",cola,((1LL<<(64-NUMBITSPREBUSQUEDAS))-1LL));
  
  for(msb=0;msb<NPREBUSQUEDAS;msb++){
    temporal=msb<<(64-NUMBITSPREBUSQUEDAS);
    for(i=0;i<256;i++){
      if(temporal<cumulativa[i])
	break;
    }
    if(i>255)
      i=255;
    tabla_min[msb]=(uint32_t) i;

    temporal|=cola;
    for(i=0;i<256;i++){
      if(temporal<cumulativa[i])
	break;
    }
    if(i>255)
      i=255;
    tabla_max[msb]=(uint32_t) i;
    if(tabla_min[msb]>tabla_max[msb])
      print_and_exit("be=%.14g, inversion en msb=%u\n",beta,msb);
  }

  //Ahora empaquetamos
  for(i=0;i<NPREBUSQUEDAS;i+=8){
    for(j=0;j<4;j++){
      paquete.vec[j]=tabla_min[i+(j<<1)]|(tabla_max[i+(j<<1)]<<8)|
	(tabla_min[i+(j<<1)+1]<<16)|(tabla_max[i+(j<<1)+1]<<24);
    }
    mi_lut[0].prebusqueda[i>>3]=paquete.palabra;    
  }

  //Comprobamos el empaquetamiento:
  for(i=0;i<NPREBUSQUEDAS;i++){
    paquete.palabra=mi_lut[0].prebusqueda[i>>3];
    j=(i&1)<<4;
    if(((paquete.vec[(i&7)>>1]>>j)&255)!=tabla_min[i])
      print_and_exit("beta=%.14g mal min de msb=%d\n",i);
    if(((paquete.vec[(i&7)>>1]>>(j+8))&255)!=tabla_max[i])
      print_and_exit("beta=%.14g mal max de msb=%d\n",i);
  }      

}

#undef TWOBRQ

void print_and_exit(const char *format, ...)
{
  va_list list;
    
  va_start(list,format);
  vfprintf(stderr,format,list);
  va_end(list);
  exit(1);
}

void escribe_luts(void)  
{
  FILE *Fout;
  char name[1024];
  sprintf(name,"LUT_for_PRNG_nbits%02d_NB%d.bin",NUMBITSPREBUSQUEDAS,NUMBETAS);
  if(NULL==(Fout=fopen(name,"wb")))
    print_and_exit("Problemas abriendo %s\n",name);

  
  if(1!=fwrite(betas,sizeof(betas),(size_t) 1,Fout))
    print_and_exit("problemas escribiendo betas\n");
  if(1!=fwrite(lut_heat_bath,sizeof(lut_heat_bath),(size_t) 1,Fout))
    print_and_exit("problemas escribiendo lut_heat_bath\n");
 
  fclose(Fout);
}

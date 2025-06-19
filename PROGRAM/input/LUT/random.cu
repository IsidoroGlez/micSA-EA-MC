#include <stdio.h>
#include "random.h"


//Las variables que definen a un generador de Luescher de 24 bits
//(solo se usan para inicializacion):
uint32_t lusch1[256];
char bit_luscher;
unsigned char l_ip,l_ip1,l_ip2;


uint64_t comprueba_semilla(uint64_t semilla)
{
  if(!semilla)
    semilla=lee_urandom();
  return semilla;
}

uint64_t lee_urandom(void)
{
  FILE *frandom;
  uint64_t output;
  frandom=fopen("/dev/urandom","r");
  fread(&output,sizeof(output),(size_t) 1,frandom);
  fclose(frandom);
  return output;
}
void Init_luxury(uint64_t semilla)
{
  int i,j;
  uint64_t temporal;
  uint32_t bobo;

  l_ip=128;    
  l_ip1=l_ip-10;    
  l_ip2=l_ip-24;  
  
  for(i=0;i<1111;i++)//Por si la semilla fuera chica...
    semilla=semilla*3202034522624059733LLU+1;

  for(j=l_ip2; j<l_ip;  j++){
    temporal=0ULL;
    for(i=0;i<3;i++){
      semilla=semilla*3202034522624059733LLU+1;
      temporal|=(semilla>>56)<<(8*i);
    }
    lusch1[j] = (uint32_t) temporal;
  }
  semilla=semilla*3202034522624059733LLU+1;
  bit_luscher=((semilla>>63)>0)?1:0;

  for(i=0;i<1111;i++)//Que corra el luxury...
    bobo=luxury();

  j=bobo;//To avoid a warning...
}

uint32_t luxury()
{
  int diferencia;
  int i;
  uint32_t output;


  for(i=0;i<389;i++){ //Esta es la luxuria
    diferencia=lusch1[l_ip1]-lusch1[l_ip2]-bit_luscher;
    lusch1[l_ip]=(diferencia+TWO24)&TWO24M1;

    bit_luscher=(diferencia&(1UL<<31))>>31; //Si diferencia<0, MSB=1

    l_ip++; l_ip1++; l_ip2++;
  }
  output=lusch1[l_ip];

  return output;
}

uint64_t Luescher_PRNG(void)
{//8 bytes de 8 llamadas consecutivas a luxury.
  
  int j;
  uint64_t output;
  uint64_t mi_byte;

  output=0;
  for(j=0;j<8;j++){
    mi_byte=(uint64_t) luxury();
    output|=(mi_byte>>16)<<8*j;
  }

  return output;
}


void Inicia_generadores_CPU(s_aleatorio_HQ_64bits * PRC, s_xoshiro256pp * xoshiro,
		       uint32_t * semilla_key, uint64_t semilla)

{
  s_aleatorio_HQ_64bits salida;
  uint32_t bobo,bobo2;
  uint64_t estado[4];
  uint64_t lectura;
  int i,j,k;
  
  semilla=comprueba_semilla(semilla);
  Init_luxury(semilla);
  
  for(i=0;i<1111;i++)//Que corra la luxuria
    bobo=luxury();

  //Parisi-Rapuano+Congruencial
  k=bobo;//Useless, just to avoid a warning.
  salida.zseed=Luescher_PRNG();

  for(k=128-61;k<128;k++)
    salida.ira[k]=Luescher_PRNG();
  
  salida.ip=128;    
  salida.ip1=salida.ip-24;    
  salida.ip2=salida.ip-55;    
  salida.ip3=salida.ip-61;
  
  for(i=0;i<2222;i++){
    _actualiza_aleatorio_HQ_escalar(salida);
  }

  PRC[0]=salida;
  
  // Ahora el xoshiro256++

  estado[0]=0ULL;
  estado[1]=0ULL;
  estado[2]=0ULL;
  estado[3]=0ULL;

  //No es aceptable una semilla de 256 bit nulos.
  while(!(estado[0]|estado[1]|estado[2]|estado[3])){
    for(j=0;j<4;j++){
      estado[j]=0ULL;
      for(i=0;i<64;i++){
        lectura=(uint64_t) luxury();
        estado[j]|=((lectura>>23)&1ULL)<<i;
      }
    }
  }

  for(j=0;j<4;j++)
    xoshiro[0].s[j]=estado[j];

  bobo2=0;
  for(i=0;i<64;i++){
    bobo=luxury();
    bobo2|=((bobo>>23)&1U)<<i;
  }
  semilla_key[0]=bobo2;
    
}

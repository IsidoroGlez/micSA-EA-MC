#include "bits-random.h"

//Init

void Init()
{
  //One could say, "wait, I can just use ib". Don't. Errors will happen for sure in other places of the program
  //if you are not always careful about what are you looking at, "clones" or "betas".
  
  int ib,ic;
  
  for(ib=0;ib<NUMBETAS;ib++)
    cual_beta_este_clon[ib]=ib;
  for(ic=0;ic<NUMBETAS;ic++)
      cual_clon_esta_beta[ic]=ic;
  acera_acumuladores();
}

void acera_acumuladores()
{
  int ic,site,tid;
  tid=0;
  for(ic=0;ic<NUMBETAS;ic++)
    for(site=0;site<MSC_V;site++){
       suma_global[tid].x= suma_global[tid].y= suma_global[tid].z= suma_global[tid].w=0;
       tid++;
    }
}

#include "../include/header.h"

void Measures(int nbits)
{
  int ibit, ibeta, ir, iclon;
  double aux;
  
  // KERNEL ENERGY
  Measure_Energy_Spins(nbits);
  Measure_Energy_Clerics(nbits);
      
  // AVERAGE ENERGY
  for(ibit=0;ibit<nbits;ibit++){
    memset(average_Energy[ibit], 0, sizeof(double)*NBETAS);   //WARNING: instantaneous measure!!
    memset(average_Energy_C[ibit], 0, sizeof(double)*NBETAS);
    
    memset(std_Energy[ibit], 0, sizeof(double)*NBETAS);   //WARNING: instantaneous measure!!
    memset(std_Energy_C[ibit], 0, sizeof(double)*NBETAS);
    for(ir=0;ir<NR;ir++){
      for(ibeta=0;ibeta<NBETAS;ibeta++){
	iclon = which_clon_this_beta[ibit][ir][ibeta];
	aux = ((double) (-h_Ener[iclon+ir*NBETAS+ibit*NRNBETAS])) / ((double) V);
	average_Energy[ibit][ibeta] += aux;
	std_Energy[ibit][ibeta] += aux*aux;

	aux = ((double) ener_C[ibit][ibeta+ir*NBETAS]) / ( ((double) V)*data.nc );
	average_Energy_C[ibit][ibeta] += aux;
	std_Energy_C[ibit][ibeta] += aux*aux;
      }
    }
  }
  
  // WRITE ENERGY MEASURES
  for(ibit=0;ibit<nbits;ibit++){
    for(ibeta=0;ibeta<NBETAS;ibeta++){
      average_Energy[ibit][ibeta] /= (double) NR;
      average_Energy_C[ibit][ibeta] /= (double) NR;

      std_Energy[ibit][ibeta] /= (double) NR;
      std_Energy_C[ibit][ibeta] /= (double) NR;
      
      std_Energy[ibit][ibeta] -= average_Energy[ibit][ibeta]*average_Energy[ibit][ibeta];
      std_Energy_C[ibit][ibeta] -= average_Energy_C[ibit][ibeta]*average_Energy_C[ibit][ibeta];

      std_Energy[ibit][ibeta] /= (double) NR-1;
      std_Energy_C[ibit][ibeta] /= (double) NR-1;

      std_Energy[ibit][ibeta] = sqrt(std_Energy[ibit][ibeta]);
      std_Energy_C[ibit][ibeta] = sqrt(std_Energy_C[ibit][ibeta]);
    }
  }

}

void Measure_Energy_Clerics(int nbits){

  unsigned long long offset = HALF_MSC_VNRNBETAS*data.nc;
  MYWORD *aguja0[2], *aguja1[2];

  int ibit, ic, ibeta, site, ir;
  
  for(ibit=0;ibit<nbits;ibit++){
    memset(ener_C[ibit], 0, sizeof(double)*NRNBETAS);
    aguja0[0] = write_CPUwalker0[ibit];
    aguja0[1] = write_CPUwalker0[ibit]+offset;
    aguja1[0] = write_CPUwalker1[ibit];
    aguja1[1] = write_CPUwalker1[ibit]+offset;
    for(ic=0;ic<data.nc;ic++)
      for(ibeta=0;ibeta<NBETAS;ibeta++)
	for(site=0;site<HALF_MSC_V;site++)
	  for(ir=0;ir<NR;ir++){
	    ener_C[ibit][ibeta+ir*NBETAS] += __builtin_popcountll(*aguja0[0]);
	    ener_C[ibit][ibeta+ir*NBETAS] += __builtin_popcountll(*aguja0[1]);
	    ener_C[ibit][ibeta+ir*NBETAS] += 2*(__builtin_popcountll(*aguja1[0]));
	    ener_C[ibit][ibeta+ir*NBETAS] += 2*(__builtin_popcountll(*aguja1[1]));
	    aguja0[0]++;
	    aguja0[1]++;
	    aguja1[0]++;
	    aguja1[1]++;
	}
  }

}

int calculate_scalar_energy(int ibit, int irep,int iclon)
{
  int energia,x,y,z,site;
  int neigh_pz,neigh_px,neigh_py;
  char * aguja;
  char spin;

  energia=0;
  aguja=uu[ibit][irep][iclon];
  site=0;
  for(z=0;z<Lz;z++){
    neigh_pz=z_p[z];
    for(y=0;y<Ly;y++){
      neigh_py=y_p[y];
      for(x=0;x<Lx;x++){
	neigh_px=x_p[x];
	spin=aguja[site];

	energia+=spin^Jx[ibit][site]^aguja[site+neigh_px];
	energia+=spin^Jy[ibit][site]^aguja[site+neigh_py];
	energia+=spin^Jz[ibit][site]^aguja[site+neigh_pz];
	
	site++;
      }
    }
  }
  energia=3*V-2*energia;
  
  return energia;
}

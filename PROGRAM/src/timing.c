#include "../include/header.h"

void measure_time(int nbits, unsigned long long t)
{
    static struct timespec clock_ini={0,0},clock_fin;
    static struct timespec CPU_ini,CPU_fin;
    static double totalspins;
    double elapsed,consumed,bobo;

    clock_gettime(CLOCK_MONOTONIC,&clock_fin);
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID,&CPU_fin);

    if (clock_ini.tv_sec>0) // la segunda vez que se llama
    { 
	elapsed =(clock_fin.tv_sec-clock_ini.tv_sec)*1.e9
	    +(clock_fin.tv_nsec-clock_ini.tv_nsec);
	consumed=(CPU_fin.tv_sec-CPU_ini.tv_sec)*1.e9
	    +(CPU_fin.tv_nsec-CPU_ini.tv_nsec);

	bobo = (double) (V*NR*nbits);
	bobo *= (double) (NBETAS*data.nc);
	bobo *= (double) t;
	totalspins = bobo;
	
	printf("%4.0fs (%.3f ps/s, %.0f%%): ",
	       elapsed*1e-9,
	       elapsed/totalspins*1e3,
	       consumed/elapsed*100);
    }
    
    clock_ini=clock_fin;
    CPU_ini=CPU_fin;
}

int check_time(unsigned long long max_time, unsigned long long dt0, unsigned long long dt1)
{
  static struct timespec clock_ini={0,0},clock_fin;
  static unsigned int total_time = 0;
  unsigned int elapsed, next_elapsed;
  double mean, std;
  

  clock_gettime(CLOCK_MONOTONIC,&clock_fin);
  
  if ( (clock_ini.tv_sec>0) && (max_time!=0) ){
    elapsed =(clock_fin.tv_sec-clock_ini.tv_sec);
    
    total_time += elapsed;

    next_elapsed = (unsigned int) elapsed*(dt1/dt0)*1.1;
    if( (((double) max_time)-total_time-next_elapsed) <= 0 )
      return 1;
  }
  
  clock_ini=clock_fin;
  
  return 0;
}

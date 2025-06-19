#define MAIN
#include <pthread.h>
#include "../include/cudamacro.h"
#include "../include/header.h"

#include "../include/fill_bits.h"

cudaStream_t *stream;

int main(int argc, char **argv){

  if ( (argc > 1) && ( (strcmp(argv[1], "-h") == 0)||(strcmp(argv[1], "--help") == 0) ) ) {
    print_help(argv[0]);
  }

  // info sample and replica
  int sample, nbits;
  
  // names of input files:
  char name_input[MAXSTR];
  char name_betas[MAXSTR];
  char name_lut[MAXSTR];
  char name_list[MAXSTR] = {'\0'};
  char name_times[MAXSTR];
  
  unsigned long long iter;
  int ibit, i;
  
  uint32_t mi_jobID;
  int device,count_devices;
  unsigned long long max_time = 0;
  
  switch(argc){
  case 10:
    snprintf(name_list,MAXSTR,"%s",argv[9]);
  case 9:
    sscanf(argv[8],"%llu",&max_time);
  case 8:
    sscanf(argv[7],"%d",&device);
    snprintf(name_times,MAXSTR,"%s",argv[6]);
    snprintf(name_lut,MAXSTR,"%s",argv[5]);
    snprintf(name_input,MAXSTR,"%s",argv[4]);
    snprintf(name_betas,MAXSTR,"%s",argv[3]);
    sscanf(argv[2],"%d",&nbits);
    sscanf(argv[1],"%d",&sample);
    break;
  default:
    print_and_exit("Usage: %s isample nbits beta.dat input.in LUT list_times device [max_time [list_samples]]\n",
		   argv[0]);
  }

  MY_CUDA_CHECK(cudaGetDeviceCount( &count_devices) );
  if((device<0)||(device>=count_devices))
    print_and_exit("device=%d count_devices=%d\n",device,count_devices);
  
  printf("count device = %d; device = %d\n", count_devices, device);
  fflush(stdout);
  
  MY_CUDA_CHECK( cudaSetDevice(device) );
  
  measure_time(nbits,0);
  check_time(max_time,1,0);
  
  if((nbits<0) | (nbits>MAXNSAMPLES))
    print_and_exit("Error: nbits out of range, 0 < bit < %d\n", MAXNSAMPLES);
  
  if(sample<0)
    print_and_exit("Error: sample out of range, 0 <= sample\n");
  
  mi_jobID = (uint32_t) sample;
  if(mi_jobID>((1<<NUMBITS_PER_JOBID)-1))
    print_and_exit("Our Counter-based PRNG  cannot accomodate more than %d samples\n",
		   (1<<NUMBITS_PER_JOBID));

  // Read INPUT DATA, BETAS and LUT 
  read_input(name_input);
  read_betas(name_betas);
  read_lut(name_lut, name_betas);
  read_list_times(name_times);
  read_list_samples(name_list, nbits);

  // Creating output path
  create_sample_path(sample, nbits);

  // Checking flag for backup and getting seeds
  check_and_prepare_simulation_backup(sample, nbits);

  printf("Betas and LUT have been read correctly.\n");
  printf("Initializating J's and Spins ...\n");
  fflush(stdout);
  
  handle_sigterm();
  create_running();

  //INIT EVERYTHING
  Init(nbits);

  //CHECK BACKUP OR NOT
  if(data.flag==0){
    backup(nbits);
  }
  
  if(write_seeds){
    write_seeds_in_file(nbits);
  }
  
  stop(&data.maxtime); // Last chance for changing nbin

  //CUDA KERNELS & STREAMS
  if(NULL==(stream = (cudaStream_t*) malloc(nbits * sizeof(cudaStream_t))))
    print_and_exit("Problems allocating streams\n");
  
  for(ibit=0;ibit<nbits;ibit++){
    MY_CUDA_CHECK( cudaStreamCreate(&(stream[ibit])) );
  }

  //COPY EVERITHING TO DEVICE MEMORY
  CopyConstOnDevice();
  CopyDataOnDevice(nbits);

  //MAIN LOOP
  send_PT_state_to_GPU(nbits); //WARNING: remanente de tener varias temperaturas
  
  iter = 0;
  // THERMALIZATION
  while(iter<data.nterm){
    iter = Monte_Carlo(mi_jobID, iter, nbits, 0);
  }

  refresh_walkers(nbits);
  
  iter = data.now;
  i = 0;
  while(data.now>=list_times[i])
    i++;
  
  while(iter < data.maxtime){

    //SYNCHRONIZE
    if( !(iter&(iter-1))){
      refresh_walkers(nbits);
    }

    MY_CUDA_CHECK( cudaDeviceSynchronize() );

    iter = Monte_Carlo(mi_jobID, iter, nbits, 1);

    //SYNCHRONIZE
    MY_CUDA_CHECK( cudaDeviceSynchronize() );

    data.now=iter;

    if(iter==list_times[i]){

      // MEASURE PROGRAM SPEED
      measure_time(nbits,(i>0)?(iter-list_times[i-1]):iter);
      
      ignore_sigterm(); 

      // Write MEASURES, CONF AND WALKERS
      copy_from_GPU_for_write(nbits);

      Measures(nbits);
      write_measures(iter, nbits);

      write_walkers(data.now, nbits, (i>0)?list_times[i-1]:0 );
      write_conf(data.now, nbits);

      stop(&data.maxtime); // para o alarga el loop
      handle_sigterm();// ya estan todos los ficheros escritos

      // CHECK IF THERE IS ENOUGH TIME
      if( (check_time(max_time,
		    (i>0)?(iter-list_times[i-1]):iter,
		    ((i+1)<data.ntimes)?list_times[i+1]-iter:0)) ||
	  ((++i)==data.ntimes) )
	break;

    }
  } //while
  
  // END SIMULATION
  ignore_sigterm(); 
  delete_running();
  
  printf("Simulation I=%d, nbits=%d finished at iteration %llu\n",
	 sample, nbits,iter);
  fflush(stdout);

  // FREE EVERYTHING
  FreeDevice(nbits);
  FreeHost(nbits);
  
  return 0;
}

int Monte_Carlo(uint32_t mi_jobID, unsigned long long t, int nbits, int flag)
{
  int ic, ir;
  unsigned long long iter;
  uint64_t parity, entropy, temporal;

  dim3 dimGrid(numBlocks,nbits,1);
  dim3 dimBlock(numThreadsPerBlock,1,1);

  if(!flag){
    iter = 2*t;
    //KERNEL MC EVEN
    parity=0;
    _fill_bits;
    d_gen_rndbits_cuda<<<numBlocks,numThreadsPerBlock>>>(HALF_MSC_V,
							 NR,
							 NBETAS,
							 s_time_and_entropy,
							 seed_keys,
							 dev_lut_heat_bath,
							 rand_wht_h);

    d_oneMCstepBN_multisample<<<dimGrid, dimBlock>>>(d_spin,
						     parity,
						     d_J,
						     d_neig,
						     HALF_MSC_V,
						     NR,
						     NBETAS,
						     d_whichclone,
						     d_bianchi_rotate,
						     rand_wht_h,
						     d_walker0, d_walker1, 0);
    iter++;
    //KERNEL MC ODD
    parity=1;
    _fill_bits;
    d_gen_rndbits_cuda<<<numBlocks,numThreadsPerBlock>>>(HALF_MSC_V,
							 NR,
							 NBETAS,
							 s_time_and_entropy,
							 seed_keys,
							 dev_lut_heat_bath,
							 rand_blk_h);
    
    d_oneMCstepBN_multisample<<<dimGrid, dimBlock>>>(d_spin,
						     parity,
						     d_J,
						     d_neig,
						     HALF_MSC_V,
						     NR,
						     NBETAS,
						     d_whichclone,
						     d_neri_rotate,
						     rand_blk_h,
						     d_walker0, d_walker1, 0);
    
  }else{
    for(ic = 0; ic < data.nc; ic++){
      //KERNEL MC EVEN
      parity=0;
      d_oneMCstepBN_multisample<<<dimGrid, dimBlock>>>(d_spin,
						       parity,
						       d_J,
						       d_neig,
						       HALF_MSC_V,
						       NR,
						       NBETAS,
						       d_whichclone,
						       d_bianchi_rotate,
						       (uint64_t *) 0,
						       d_walker0, d_walker1, ic);
    
      //KERNEL MC ODD
      parity=1;    
      d_oneMCstepBN_multisample<<<dimGrid, dimBlock>>>(d_spin,
						       parity,
						       d_J,
						       d_neig,
						       HALF_MSC_V,
						       NR,
						       NBETAS,
						       d_whichclone,
						       d_neri_rotate,
						       (uint64_t *) 0,
						       d_walker0, d_walker1, ic);
      
    }
#ifndef NONWALK
    parity = 0;
    d_Walk<<<dimGrid, dimBlock>>>(parity, d_neig,
					  HALF_MSC_V,
					  NR,
					  NBETAS,
					  data.nc,
					  d_bianchi_rotate,
					  d_walker0,
					  d_walker1);
    
    parity = 1;
    d_Walk<<<dimGrid, dimBlock>>>(parity, d_neig,
					  HALF_MSC_V,
					  NR,
					  NBETAS,
					  data.nc,
					  d_neri_rotate,
					  d_walker0,
					  d_walker1);
#endif
    
  }

  t++;
  
  return t;
  
}

void Measure_Energy_Spins (int nbits)
{
  int ibit;
  
  MY_CUDA_CHECK( cudaMemset( d_Ener, 0, nbits*NRNBETAS*sizeof(int) ) );

  for(ibit=0;ibit<nbits;ibit++){
    d_computeEne<<<numBlocks,numThreadsPerBlock,numThreadsPerBlock*sizeof(MYWORD),stream[ibit]>>>(h_spin[ibit][0],
      h_spin[ibit][1],
      d_Ener+(ibit*NRNBETAS),
      d_sum, 0,
      d_J[ibit],
      d_neig[ibit],
      HALF_MSC_V,
      NR, NBETAS,
      d_bianchi_index,
      d_neri_index);
  }

  MY_CUDA_CHECK( cudaMemcpy(h_Ener, d_Ener, sizeof(int)*NRNBETAS*nbits,
			    cudaMemcpyDeviceToHost) );
}

// DEVICE MEMORY

void refresh_walkers(int nbits){
  static int first = 1;

  int rc;
  static pthread_t ptidic[MAXNSAMPLES];
  
  if(!(first*data.flag)){
    for(int g=0; g<nbits; g++) {
      
      if(!first){
	if((rc=pthread_join(ptidic[g],NULL))<0) {
	  print_and_exit("Error in join thread %d:\n",rc);
	}
      }
      
      MY_CUDA_CHECK( cudaMemcpy(h_walker0[g][0],CPUwalker0[g],
				sizeof(MYWORD)*HALF_MSC_VNRNBETAS*data.nc,cudaMemcpyHostToDevice) );
      MY_CUDA_CHECK( cudaMemcpy(h_walker0[g][1],CPUwalker0[g]+HALF_MSC_VNRNBETAS*data.nc,
				sizeof(MYWORD)*HALF_MSC_VNRNBETAS*data.nc,cudaMemcpyHostToDevice) );
      MY_CUDA_CHECK( cudaMemcpy(h_walker1[g][0],CPUwalker1[g],
				sizeof(MYWORD)*HALF_MSC_VNRNBETAS*data.nc,cudaMemcpyHostToDevice));
      MY_CUDA_CHECK( cudaMemcpy(h_walker1[g][1],CPUwalker1[g]+HALF_MSC_VNRNBETAS*data.nc,
				sizeof(MYWORD)*HALF_MSC_VNRNBETAS*data.nc,cudaMemcpyHostToDevice) );

      memcpy(write_random_PRC_C[g], random_PRC_C[g],sizeof(s_aleatorio_HQ_64bits)*NR);
      memcpy(write_random_xoshiro256pp_C[g], random_xoshiro256pp_C[g], sizeof(s_xoshiro256pp)*NR);
    }
  }

  first = 0;
  
  for(int g=0;g<nbits;g++){
    if(pthread_create(ptidic+g,NULL,InitWalkers,p2icds[g])<0) {
      print_and_exit("Error creating thread for initializing daemons/walkers!\n");
    }
  }
  
}

void CopyConstOnDevice(void)
{

  size_t total;
  int scratch;
  
  MY_CUDA_CHECK( cudaMemcpyToSymbol(d_deltaE,&off[0],sizeof(int),0,cudaMemcpyHostToDevice) );
  MY_CUDA_CHECK( cudaMemcpyToSymbol(d_deltaN,&off[1],sizeof(int),0,cudaMemcpyHostToDevice) );	
  MY_CUDA_CHECK( cudaMemcpyToSymbol(d_deltaU,&off[2],sizeof(int),0,cudaMemcpyHostToDevice) );
  scratch=-off[0];
  MY_CUDA_CHECK( cudaMemcpyToSymbol(d_deltaO,&scratch,sizeof(int),0,cudaMemcpyHostToDevice) );
  scratch=-off[1];
  MY_CUDA_CHECK( cudaMemcpyToSymbol(d_deltaS,&scratch,sizeof(int),0,cudaMemcpyHostToDevice) );	
  scratch=-off[2];
  MY_CUDA_CHECK( cudaMemcpyToSymbol(d_deltaD,&scratch,sizeof(int),0,cudaMemcpyHostToDevice) );

  total = NBETAS*sizeof(s_lut_heat_bath);
  MY_CUDA_CHECK( cudaMalloc((void **) &dev_lut_heat_bath,total));
  MY_CUDA_CHECK( cudaMemcpy( dev_lut_heat_bath, h_LUT, total , cudaMemcpyHostToDevice) );  

  //  WARNING: Por qué no se copian las J's a una constante ?
  //  WARNING: Lo mismo para los index y las rotaciones blancas y negras...
  
}

void CopyDataOnDevice(int nbits){

  size_t  total;
  int ibit;
  
  Vicini *vicini;
  MYWORD jtemp[HALF_MSC_V];
  int neigtemp[HALF_MSC_V];
  int d, i, j;
  
  total = sizeof(unsigned int)*HALF_MSC_V;

  MY_CUDA_CHECK( cudaMalloc(&d_bianchi_index, total) );
  MY_CUDA_CHECK( cudaMalloc(&d_neri_index, total) );	
  MY_CUDA_CHECK( cudaMemcpy(d_bianchi_index, bianchi_index, total,cudaMemcpyHostToDevice) );
  MY_CUDA_CHECK( cudaMemcpy(d_neri_index, neri_index, total, cudaMemcpyHostToDevice) );

  total = sizeof(unsigned char)*HALF_MSC_V;
  
  MY_CUDA_CHECK( cudaMalloc(&d_bianchi_rotate, total) );
  MY_CUDA_CHECK( cudaMalloc(&d_neri_rotate, total) );
  MY_CUDA_CHECK( cudaMemcpy(d_bianchi_rotate, bianchi_rotate, total,cudaMemcpyHostToDevice) );
  MY_CUDA_CHECK( cudaMemcpy(d_neri_rotate, neri_rotate, total, cudaMemcpyHostToDevice) );

  for(ibit=0;ibit<nbits;ibit++){
    h_spin[ibit] = (MYWORD **)mmcuda((void ***)&d_spin[ibit],2,HALF_MSC_VNRNBETAS,sizeof(MYWORD),1);

    MY_CUDA_CHECK( cudaMemcpy(h_spin[ibit][0], h_MSC_u_even[ibit],
			      sizeof(MYWORD)*HALF_MSC_VNRNBETAS, cudaMemcpyHostToDevice) );
    MY_CUDA_CHECK( cudaMemcpy(h_spin[ibit][1], h_MSC_u_odd[ibit],
			      sizeof(MYWORD)*HALF_MSC_VNRNBETAS, cudaMemcpyHostToDevice) );
  
    h_J[ibit] = (MYWORD **)mmcuda((void ***)&d_J[ibit],2*DEGREE,HALF_MSC_V,sizeof(MYWORD),1);
    h_neig[ibit] = (int **)mmcuda((void ***)&d_neig[ibit],2*DEGREE,HALF_MSC_V,sizeof(int),1);


    h_walker0[ibit]=(MYWORD **)mmcuda((void ***)&d_walker0[ibit],2,HALF_MSC_VNRNBETAS*data.nc,sizeof(MYWORD),1);
    h_walker1[ibit]=(MYWORD **)mmcuda((void ***)&d_walker1[ibit],2,HALF_MSC_VNRNBETAS*data.nc,sizeof(MYWORD),1);   
    
    for(d=0; d<2; d++) {
      vicini = (d==0)?viciniB:viciniN;
      for(i=0; i<DEGREE; i++) {
      for(j=0; j<HALF_MSC_V; j++) {
	jtemp[j]=vicini[j].J[ibit][i];
	neigtemp[j]=vicini[j].neig[i];
      }
      MY_CUDA_CHECK( cudaMemcpy(h_J[ibit][i+d*DEGREE], jtemp,
				sizeof(MYWORD)*HALF_MSC_V, cudaMemcpyHostToDevice) );    
      MY_CUDA_CHECK( cudaMemcpy(h_neig[ibit][i+d*DEGREE], neigtemp,
				sizeof(int)*HALF_MSC_V, cudaMemcpyHostToDevice) );
      }
    }
  
    MY_CUDA_CHECK( cudaMalloc(&d_whichclone[ibit],sizeof(unsigned int)*NRNBETAS) );
  }

  if(NULL==(h_Ener = (int *) malloc(sizeof(int)*NRNBETAS*nbits)))
    print_and_exit("Problems allocating h_Ener\n");
  MY_CUDA_CHECK( cudaMalloc(&d_Ener,sizeof(int)*NRNBETAS*nbits) );
  
  MY_CUDA_CHECK( cudaMalloc(&d_sum,sizeof(int)*MAX) );
  MY_CUDA_CHECK( cudaMemcpy(d_sum, sum, sizeof(sum), cudaMemcpyHostToDevice) );

  total = 3*HALF_MSC_VNRNBETAS*sizeof(uint64_t);
  MY_CUDA_CHECK(cudaMalloc(&rand_wht_h, total));
  MY_CUDA_CHECK(cudaMalloc(&rand_blk_h, total));
  
}

void inline copy_from_GPU_for_write(int nbits)
{
  for(int ibit=0;ibit<nbits;ibit++){
    MY_CUDA_CHECK( cudaMemcpy(h_MSC_u_even[ibit], h_spin[ibit][0],
			      sizeof(MYWORD)*HALF_MSC_VNRNBETAS, cudaMemcpyDeviceToHost) );
    MY_CUDA_CHECK( cudaMemcpy(h_MSC_u_odd[ibit], h_spin[ibit][1],
			      sizeof(MYWORD)*HALF_MSC_VNRNBETAS, cudaMemcpyDeviceToHost) );
	
    MY_CUDA_CHECK( cudaMemcpy(write_CPUwalker0[ibit], h_walker0[ibit][0],
			      sizeof(MYWORD)*HALF_MSC_VNRNBETAS*data.nc,cudaMemcpyDeviceToHost) );
    MY_CUDA_CHECK( cudaMemcpy(write_CPUwalker0[ibit]+HALF_MSC_VNRNBETAS*data.nc, h_walker0[ibit][1],
			      sizeof(MYWORD)*HALF_MSC_VNRNBETAS*data.nc,cudaMemcpyDeviceToHost) );
    MY_CUDA_CHECK( cudaMemcpy(write_CPUwalker1[ibit], h_walker1[ibit][0],
			      sizeof(MYWORD)*HALF_MSC_VNRNBETAS*data.nc,cudaMemcpyDeviceToHost));
    MY_CUDA_CHECK( cudaMemcpy(write_CPUwalker1[ibit]+HALF_MSC_VNRNBETAS*data.nc, h_walker1[ibit][1],
			      sizeof(MYWORD)*HALF_MSC_VNRNBETAS*data.nc,cudaMemcpyDeviceToHost) );
  }
}

void send_PT_state_to_GPU(int nbits)
{

  unsigned int h_whichclone[NR*NBETAS];
  int r, ib, ibit;

  for(ibit=0;ibit<nbits;ibit++){
    for(r=0;r<NR;r++)
      for(ib=0;ib<NBETAS;ib++)
	h_whichclone[r+NR*ib] = (unsigned int) which_beta_this_clon[ibit][r][ib];
  
    MY_CUDA_CHECK( cudaMemcpy(d_whichclone[ibit], h_whichclone,
			      sizeof(unsigned int)*NRNBETAS, cudaMemcpyHostToDevice) );
  }
}

void FreeDevice(int nbits){

  MY_CUDA_CHECK(cudaFree(d_bianchi_index));
  MY_CUDA_CHECK(cudaFree(d_neri_index));

  MY_CUDA_CHECK(cudaFree(d_bianchi_rotate));
  MY_CUDA_CHECK(cudaFree(d_neri_rotate));

  MY_CUDA_CHECK(cudaFree(d_sum));

  MY_CUDA_CHECK(cudaFree(dev_lut_heat_bath));

  MY_CUDA_CHECK(cudaFree(rand_wht_h));
  MY_CUDA_CHECK(cudaFree(rand_blk_h));

  MY_CUDA_CHECK(cudaFree(d_Ener));
  
  for(int ibit=0;ibit<nbits;ibit++){
    MY_CUDA_CHECK(cudaFree(d_whichclone[ibit]));
  }
}

void FreeHost(int nbits)
{
  for(int g=0;g<nbits;g++){
    free(CPUwalker0[g]);
    free(CPUwalker1[g]);
    free(write_CPUwalker0[g]);
    free(write_CPUwalker1[g]);    
  }

  free(stream);
  free(h_Ener);
  free(list_times);
 
}


// KERNEL RNG
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

#define _update_key {\
    key[0]+=PHILOX_W32_0; \
    key[1]+=PHILOX_W32_1;}

#define _update_state {\
  lo0=PHILOX_M4x32_0*v[0];                      \
  hi0=__umulhi(PHILOX_M4x32_0,v[0]);            \
  lo1=PHILOX_M4x32_1*v[2];                      \
  hi1=__umulhi(PHILOX_M4x32_1,v[2]);            \
  v[0]=hi1^v[1]^key[0];                         \
  v[1]=lo1;                                     \
  v[2]=hi0^v[3]^key[1];                         \
  v[3]=lo0;}

#define _philox_4x32_10 {\
  _update_state; \
  _update_key;   \
  _update_state; \
  _update_key;   \
  _update_state; \
  _update_key;   \
  _update_state; \
  _update_key;   \
  _update_state; \
  _update_key;   \
  _update_state; \
  _update_key;   \
  _update_state; \
  _update_key;   \
  _update_state; \
  _update_key;   \
  _update_state; \
  _update_key;   \
  _update_state;}

#define _obten_aleatorio_contador {					\
    uint32_t hi0,hi1,lo0,lo1;						\
    uint32_t key[2];							\
    key[0]=tid; key[1]=useed;	\
    v[0]=time_and_entropy.x;						\
    v[1]=time_and_entropy.y;						\
    v[2]=time_and_entropy.z;						\
    v[3]=time_and_entropy.w;						\
    _philox_4x32_10;							\
    time_and_entropy.x++; /* contador interno*/				\
    hi1=v[0]^v[2];							\
    lo1=v[1]^v[3];							\
    lsb[0]=((v[1]<<17)|(v[0]>>15))^hi1^((hi1<<21)|lo1>>11);		\
    lsb[1]=((v[0]<<17)|(v[1]>>15))^lo1^(lo1<<21);			\
    lsb[2]=(hi1<<28)|(lo1>>4);						\
    lsb[3]=(lo1<<28)|(hi1>>4);}


#define _bisection(Rmsb,Rlsb) {						\
    const uint32_t msb=Rmsb>>(32-NUMBITSPREBUSQUEDAS);			\
    uint32_t min,max,decision;						\
    decision=(uint32_t) (~(((Rmsb==dev_umbrales[254].x)?(Rlsb<dev_umbrales[254].y):(Rmsb<dev_umbrales[254].x))-1)); \
    min=dev_prebusqueda[msb]&255;					\
    max=dev_prebusqueda[msb]>>8;					\
    min=(decision&min)|((~decision)&255);				\
    max=(decision&max)|((~decision)&255);				\
    while((max-min)>1){							\
      npr=(max+min)>>1;							\
      decision=(uint32_t) (~(((Rmsb==dev_umbrales[npr].x)?(Rlsb<dev_umbrales[npr].y):(Rmsb<dev_umbrales[npr].x))-1)); \
      max=(decision&npr)|((~decision)&max);				\
      min=(decision&min)|((~decision)&npr);				\
    }									\
  decision=(uint32_t) (~(((Rmsb==dev_umbrales[min].x)?(Rlsb<dev_umbrales[min].y):(Rmsb<dev_umbrales[min].x))-1)); \
  npr=(decision&min)|(~decision&max);}

#define _get_rnd_word(bb) {				\
    uint32_t v[4],lsb[4];				\
    _obten_aleatorio_contador;				\
    const uint32_t rot=(lsb[0]+v[1]+lsb[2]+v[3])>>27;	\
    _bisection(v[0],lsb[0]);				\
    bb=npr;						\
    _bisection(v[1],lsb[1]);				\
    bb|=npr<<8;						\
    _bisection(v[2],lsb[2]);				\
    bb|=npr<<16;						\
    _bisection(v[3],lsb[3]);				\
    bb|=npr<<24;						\
    bb=(bb<<rot)|(bb>>(32-rot));}


__global__ void d_gen_rndbits_cuda(int n,
				   int nrep,
				   int nclone,
				   s_time s_time_and_entropy,
				   s_keys s_key,
				   s_lut_heat_bath *dev_lut_heat_bath,
				   uint64_t *rand_d) {

  const unsigned int tid = threadIdx.x+blockDim.x*blockIdx.x;
  const unsigned int tthreads = gridDim.x*blockDim.x;
  
  MYWORD scra;
  MYWORD b0,b1,b2;
  uint32_t b32t;

  int npr;

  unsigned int ibeta = tid / (n*nrep); //whichclone[iclone*nrep+whichrep];
  int ir=(tid-(ibeta*n*nrep))%nrep;
  
  //Select my counter and philox parameter
  uint4 time_and_entropy = s_time_and_entropy.vec[ir];
  uint32_t useed = s_key.my_key[ir];
  
  if (ibeta >= nclone) {
    return;
  }

  time_and_entropy.z += tid; 
  const uint2 *__restrict__ dev_umbrales = reinterpret_cast<const uint2 *>(dev_lut_heat_bath[ibeta].umbrales);
  const unsigned short *__restrict__ dev_prebusqueda = reinterpret_cast<const unsigned short *>(dev_lut_heat_bath[ibeta].prebusqueda);

  rand_d += tid;

  _get_rnd_word(b32t);
  scra=b32t;
  b0=scra<<32;
  _get_rnd_word(b32t);
  b0|=b32t;
  _get_rnd_word(b32t);
  scra=b32t;
  b1=scra<<32;
  _get_rnd_word(b32t);
  b1|=b32t;
  _get_rnd_word(b32t);
  scra=b32t;
  b2=scra<<32;
  _get_rnd_word(b32t);
  b2|=b32t;
  
  rand_d[0*tthreads] = b0;
  rand_d[1*tthreads] = b1;
  rand_d[2*tthreads] = b2;

  return;
}

// KERNEL ENERGY

#define MASK_E (0x1111111111111111ull)
#define MASK_N (0x000F000F000F000Full)
#define MASK_U (0x000000000000FFFFull)
#define MASK_O (0x8888888888888888ull)
#define MASK_S (0xF000F000F000F000ull)
#define MASK_D (0xFFFF000000000000ull)

#define NMASK_E (~MASK_E)
#define NMASK_N (~MASK_N)
#define NMASK_U (~MASK_U)
#define NMASK_O (~MASK_O)
#define NMASK_S (~MASK_S)
#define NMASK_D (~MASK_D)

__device__ MYWORD RotateE(MYWORD op) {
  return (op & MASK_E) << 3 | (op & NMASK_E) >> 1;
}

__device__ MYWORD RotateN(MYWORD op) {
  return (op & MASK_N) << 12 | (op & NMASK_N) >> 4;
}

__device__ MYWORD RotateU(MYWORD op) {
  return (op & MASK_U) << 48 | (op & NMASK_U) >> 16;
}

__device__ MYWORD RotateO(MYWORD op) {
  return (op & MASK_O) >> 3 | (op & NMASK_O) << 1;
}

__device__ MYWORD RotateS(MYWORD op) {
  return (op & MASK_S) >> 12 | (op & NMASK_S) << 4;
}

__device__ MYWORD RotateD(MYWORD op) {
  return (op & MASK_D) >> 48 | (op & NMASK_D) << 16;
}


__global__ void
__launch_bounds__(1024, 1)
  d_computeEne(MYWORD * __restrict__ newSpin,
	       MYWORD * __restrict__ oldSpin,
	       int * __restrict__ ene, int *sum,
	       const int dir,
	       MYWORD **J,
	       int **neig,
	       int n,
	       int nrep, int nclone,
	       unsigned int *p2p, unsigned int *p2np) {

  const unsigned int tid = threadIdx.x+blockDim.x*blockIdx.x;
  const unsigned int tthreads = gridDim.x*blockDim.x;
  const unsigned int iclone=tid/n/nrep;
  const unsigned int cloneioff=iclone*n*nrep;
  int whichrep=(tid-(iclone*n*nrep))%nrep;
  int whichsite=(tid-(iclone*n*nrep))/nrep;

  if(iclone>=nclone || whichsite>=n) return;

  int punto_n;
  int le=0;
  union Hack {unsigned long long lungo;  unsigned short corto[4];} hack;
  for(; whichsite<n; whichsite+=(tthreads/(nclone))) {
    const int punto_c=(whichsite)*nrep;
    MYWORD scra;

    punto_n=neig[DEGREE*dir+0][whichsite];
    scra=oldSpin[cloneioff+punto_n*nrep+whichrep];
    if((p2np[punto_n]-p2p[whichsite])!=d_deltaE){
      scra=RotateE(scra);
    }
    hack.lungo=newSpin[cloneioff+punto_c+whichrep]^(J[0+DEGREE*dir][whichsite]^scra);
    le+=(sum[hack.corto[0]]+sum[hack.corto[1]]+sum[hack.corto[2]]+sum[hack.corto[3]]);

    punto_n=neig[DEGREE*dir+1][whichsite];
    scra=oldSpin[cloneioff+punto_n*nrep+whichrep];
    if((p2np[punto_n]-p2p[whichsite])!=d_deltaN)scra=RotateN(scra);
    hack.lungo=newSpin[cloneioff+punto_c+whichrep]^(J[1+DEGREE*dir][whichsite]^scra);
    le+=(sum[hack.corto[0]]+sum[hack.corto[1]]+sum[hack.corto[2]]+sum[hack.corto[3]]);

    punto_n=neig[DEGREE*dir+2][whichsite];
    scra=oldSpin[cloneioff+punto_n*nrep+whichrep];
    if((p2np[punto_n]-p2p[whichsite])!=d_deltaU)scra=RotateU(scra);
    hack.lungo=newSpin[cloneioff+punto_c+whichrep]^(J[2+DEGREE*dir][whichsite]^scra);
    le+=(sum[hack.corto[0]]+sum[hack.corto[1]]+sum[hack.corto[2]]+sum[hack.corto[3]]);

    punto_n=neig[DEGREE*dir+3][whichsite];
    scra=oldSpin[cloneioff+punto_n*nrep+whichrep]^J[3+DEGREE*dir][whichsite];
    if((p2np[punto_n]-p2p[whichsite])!=d_deltaO){
      scra=RotateO(scra);
    }
    hack.lungo=scra^newSpin[cloneioff+punto_c+whichrep];
    le+=(sum[hack.corto[0]]+sum[hack.corto[1]]+sum[hack.corto[2]]+sum[hack.corto[3]]);

    punto_n=neig[DEGREE*dir+4][whichsite];
    scra=oldSpin[cloneioff+punto_n*nrep+whichrep]^J[4+DEGREE*dir][whichsite];
    if((p2np[punto_n]-p2p[whichsite])!=d_deltaS)scra=RotateS(scra);
    hack.lungo=newSpin[cloneioff+punto_c+whichrep]^scra;
    le+=(sum[hack.corto[0]]+sum[hack.corto[1]]+sum[hack.corto[2]]+sum[hack.corto[3]]);

    punto_n=neig[DEGREE*dir+5][whichsite];
    scra=oldSpin[cloneioff+punto_n*nrep+whichrep]^J[5+DEGREE*dir][whichsite];
    if((p2np[punto_n]-p2p[whichsite])!=d_deltaD)scra=RotateD(scra);
    hack.lungo=newSpin[cloneioff+punto_c+whichrep]^scra;
    le+=(sum[hack.corto[0]]+sum[hack.corto[1]]+sum[hack.corto[2]]+sum[hack.corto[3]]);
  }
  int cubetti = 6*BITSINMYWORD - 2*le;
  atomicAdd(ene + (whichrep*nclone)+iclone, cubetti);    
}

// KERNEL FOR MC

#if defined(USE_LDG)
#define LDG(x) (__ldg(&(x)))
#else
#define LDG(x) (x)
#endif

#ifdef MASSIMO
__device__ void Sum(MYWORD * sum_bit3, MYWORD *sum_bit2, MYWORD * sum_bit1,
		    MYWORD num1_bit3, MYWORD num1_bit2, MYWORD num1_bit1,
		    MYWORD num2_bit3, MYWORD num2_bit2, MYWORD num2_bit1)
{
  MYWORD carry_bit1,carry_bit2;
  sum_bit1[0] = num1_bit1 ^ num2_bit1;
  carry_bit1 = num1_bit1 & num2_bit1;
  sum_bit2[0] = num1_bit2 ^ num2_bit2 ^ carry_bit1;
  carry_bit2 = (num1_bit2 & num2_bit2) | (num1_bit2 & carry_bit1) | (num2_bit2 & carry_bit1);
  sum_bit3[0] = num1_bit3 ^ num2_bit3 ^ carry_bit2;
}
#endif
#ifdef VICTOR
#define sum_2_2(aa0,aa1,aa2,bb0,bb1,cc0,cc1){	\
    aa0 = bb0 & cc0 ;				\
    aa2 = bb1 ^ cc1 ;				\
    aa1 = aa2 ^ aa0 ;				\
    aa2 = (bb1 & cc1) ^ (aa2 & aa0) ;		\
    aa0 = bb0 ^ cc0 ;}
#endif

__global__ void d_oneMCstepBN_multisample(MYWORD **spinAll[],
					  const int dir,
					  MYWORD **JAll[],
					  int **neigAll[],
					  int n,
					  int nrep,
					  int nclone,
					  unsigned int **whichcloneAll,
					  const unsigned char *rotate,
					  uint64_t *rand_h,
					  MYWORD **walker0All[],
					  MYWORD **walker1All[], int ic) {

  const unsigned int sampleId = blockIdx.y;

  MYWORD * __restrict__ newSpin = dir ? spinAll[sampleId][1] : spinAll[sampleId][0];
  MYWORD * __restrict__ oldSpin = dir ? spinAll[sampleId][0] : spinAll[sampleId][1];

  const MYWORD **J = (const MYWORD **) JAll[sampleId] + dir*DEGREE;
  const int **neig = (const int **) neigAll[sampleId];

  MYWORD *walker0=dir ? walker0All[sampleId][1]:walker0All[sampleId][0];
  MYWORD *walker1=dir ? walker1All[sampleId][1]:walker1All[sampleId][0];

  walker0 += ic*HALF_MSC_VNRNBETAS;
  walker1 += ic*HALF_MSC_VNRNBETAS;
  
  const unsigned int *whichclone = whichcloneAll[sampleId];

  MYWORD flip;
  
  const unsigned int tid = threadIdx.x+blockDim.x*blockIdx.x;
  const unsigned int tthreads = gridDim.x*blockDim.x;

  const unsigned int iclone = tid/(n*nrep);
  const unsigned int cloneioff = iclone*n*nrep;

  int tidModNNRep = tid - (iclone*n*nrep);

  int whichrep  = tidModNNRep % nrep;
  int whichsite = tidModNNRep / nrep;

  if(iclone>=nclone || whichsite>=n) return;

  unsigned int ibeta=whichclone[iclone*nrep+whichrep];

  if(rand_h!=(uint64_t *)0) { 
    rand_h += ibeta*(n*nrep) + tidModNNRep;
  }
  
  neig += DEGREE*dir;

  oldSpin += cloneioff + whichrep;

  const MYWORD spintbu = LDG(newSpin[cloneioff + whichsite*nrep + whichrep]);

  const int punto_n0 = neig[0][whichsite];
  const int punto_n1 = neig[1][whichsite];
  const int punto_n2 = neig[2][whichsite];
  const int punto_n3 = neig[3][whichsite];
  const int punto_n4 = neig[4][whichsite];
  const int punto_n5 = neig[5][whichsite];

  MYWORD scra0 = oldSpin[punto_n0*nrep];
  MYWORD scra1 = oldSpin[punto_n1*nrep];
  MYWORD scra2 = oldSpin[punto_n2*nrep];
  MYWORD scra3 = oldSpin[punto_n3*nrep] ^ J[3][whichsite]; // negations moved into spintbu
  MYWORD scra4 = oldSpin[punto_n4*nrep] ^ J[4][whichsite];
  MYWORD scra5 = oldSpin[punto_n5*nrep] ^ J[5][whichsite];

  unsigned char rot = rotate[whichsite];
  if (rot & 0x01) scra0 = RotateE(scra0);
  if (rot & 0x02) scra1 = RotateN(scra1);
  if (rot & 0x04) scra2 = RotateU(scra2);
  if (rot & 0x08) scra3 = RotateO(scra3);
  if (rot & 0x10) scra4 = RotateS(scra4);
  if (rot & 0x20) scra5 = RotateD(scra5);
		
  MYWORD f0 = spintbu ^ J[0][whichsite] ^ scra0; // negations moved into spintbu
  MYWORD f1 = spintbu ^ J[1][whichsite] ^ scra1;
  MYWORD f2 = spintbu ^ J[2][whichsite] ^ scra2;

  //Algebra del Metropolis
  MYWORD k1 = f0 ^ f1;
  MYWORD k2 = f0 & f1;
  MYWORD j1 = k1 ^ f2;
  MYWORD k3 = k1 & f2;
  MYWORD j2 = k2 ^ k3;
		
  f0 = scra3 ^ spintbu;
  f1 = scra4 ^ spintbu;
  f2 = scra5 ^ spintbu;

  k1 = f0 ^ f1;
  k2 = f0 & f1;
  k3 = k1 & f2;
  MYWORD j3 = k1 ^ f2;
  MYWORD j4 = k2 ^ k3;
  if(rand_h!=(uint64_t *)0) { 
    MYWORD b0 = rand_h[0*tthreads];
    MYWORD b1 = rand_h[1*tthreads];
    MYWORD b2 = rand_h[2*tthreads];

    b1=b0&b1;
    b2=b1&b2;

    MYWORD id2=b1;
    MYWORD id1=(b0^b1)|b2;

    MYWORD j2ORj4=j2|j4;
    flip=(j1&j3&id1) | ((j2ORj4|id2)&(j1|j3|id1))|((j2&j4)|(j2ORj4&id2)); 

  } else {
#if defined(MASSIMO)
    MYWORD sigma0, sigma1, sigma2, ris0, ris1, OD0, OD1, D0, D1, mask;
    //G. Parisi boolean operations:
    sigma0 = j1 & j3 ;
    sigma2 = j2 ^ j4 ;
    sigma1 = sigma2 ^ sigma0 ;
    sigma2 = (j2 & j4) ^ (sigma2 & sigma0) ;
    sigma0 = j1 ^ j3 ; //so far so good
    OD1=walker1[cloneioff + whichsite*nrep + whichrep];
    OD0=walker0[cloneioff + whichsite*nrep + whichrep];
    Sum(&flip,&ris1,&ris0, ~sigma2,~sigma1,~sigma0, ~(~sigma2&(sigma0)&(sigma1)), OD1, OD0);

    D1=(flip&ris1)|((~flip)&OD1);
    D0=(flip&ris0)|((~flip)&OD0);
    /* Corrections */
    mask=(sigma0&(OD0^OD1));
    flip^=mask;
    D0^=mask;
    D1^=(sigma0&(~(OD1^sigma2)));
    /* End corrections */
    walker1[cloneioff + whichsite*nrep + whichrep]=D1;
    walker0[cloneioff + whichsite*nrep + whichrep]=D0;
    //  flip=(j1^j3)&(j2^j4); /* this was for the "pure MC" */
#elif defined(VICTOR)
    //Boolean proposed by V. Martin-Mayor
    MYWORD sigma0, sigma1, sigma2, sumando0, sumando1, sumar, cambiarR, cambiarS, carry, restando0, restando1, D0R, D1R, D0S, D1S, D0, D1, overflow, restar,neutral;
    //Walker bits taken correctly
    D0=walker0[cloneioff + whichsite*nrep + whichrep];
    D1=walker1[cloneioff + whichsite*nrep + whichrep];
    //Start boolean algebra:
    //
    // msb at right, lsb at left
    // First (j1,j2)+(j3,j4)=(sigma0,sigma1,sigma2);
    sum_2_2(sigma0,sigma1,sigma2,j1,j2,j3,j4);

    //Classification
    sumar=sigma2; //Walker get energy
    neutral=(~sigma2)&sigma1&sigma0; //Walker does not change
    restar=(~(neutral|sumar)); //walker loses energy

    //Subcases for the sum
    //Calcolate suma_global-3 (consider it is possitive,
    //because we will ise a &sumar at the end)
    //If there is overflow, the waker cannot get more energy

    sumando1=sigma1|sigma0;
    sumando0=~sigma0;
    sum_2_2(D0S,D1S,overflow,D0,D1,sumando0,sumando1);
    D0S=((D0&overflow)|(D0S&(~overflow)));
    D1S=((D1&overflow)|(D1S&(~overflow)));
    cambiarS=(~overflow);

    //Subcase for the difference
    //First check if we can accept the spin-flip
    //We can consider tht 
    //sigma0+2*sigma1 \in {0,1,2}
    cambiarR=(D1&sigma1)|((D1|sigma1)&(D0|sigma0));

    //Now we need the two's complement according to cambiarR=1;
    //sigma0+2*sigma1 =0 --->1=(1,0) (equivalent to -3 if cambiarR=1)
    //sigma0+2*sigma1 =1 --->2=(0,1) (equivalent to -2 if cambiarR=1)
    //sigma0+2*sigma1 =2 --->3=(1,1) (equivalent to -1 if cambiar R=1)
    restando1=sigma1|sigma0;
    restando0=~sigma0;

    //with the two's complement we can continue
    D0R=((restando0^D0)&cambiarR)|(D0&(~cambiarR));
    carry=restando0&D0;
    D1R=((carry^restando1^D1)&cambiarR)|(D1&(~cambiarR));

    //Chose correct subcase

    flip=neutral|(sumar&cambiarS)|(restar&cambiarR);
    D0=(neutral&D0)|(sumar&D0S)|(restar&D0R);
    D1=(neutral&D1)|(sumar&D1S)|(restar&D1R);
    walker1[cloneioff + whichsite*nrep + whichrep]=D1;
    walker0[cloneioff + whichsite*nrep + whichrep]=D0;
#elif defined(ISI)
    //Boolean algrebra proposed by I. Gonzalez-Adalid Pemartin
    MYWORD J1XORJ3, J2XORJ4, restar, sum0, sum1, D0, D1, D0_new, D1_new, carry;

    //Walker bits taken correctly:
    D0=walker0[cloneioff + whichsite*nrep + whichrep];
    D1=walker1[cloneioff + whichsite*nrep + whichrep];

    //Star boolean algrebra:
    J1XORJ3 = j1^j3;
    J2XORJ4 = j2^j4;

    // Determine if we are going to perform a sum or a difference
    restar = (J2XORJ4&( (~j1)&(~j3) )) | ( (~j2)&(~j4) );

    //Calculate the term in the operation
    sum0 = ~J1XORJ3;
    sum1 = (~J2XORJ4) & ( J1XORJ3 | ( sum0&(~(j1^j2)) ) );

    //Value of sum1 is different for the sum and the difference (apply some two's complement)
    sum1 = ( (~restar)&sum1 ) | ( restar&(sum1^sum0) );

    //We make the operation
    D0_new = D0^sum0;
    k1 = D0&sum0;
    k2 = D1^sum1;
    D1_new = k1^k2;
    k3 = D1&sum1;
    carry =  k3 | (k2&k1) ;
    carry = carry^restar;

    //Carry gives the aceptance of the flip
    flip = ~carry ;

    D0 = (flip&D0_new) | ( (~flip)&D0 ) ;
    D1 = (flip&D1_new) | ( (~flip)&D1 ) ;
    walker1[cloneioff + whichsite*nrep + whichrep]=D1;
    walker0[cloneioff + whichsite*nrep + whichrep]=D0;
#endif
  }
  
  newSpin[cloneioff + whichsite*nrep + whichrep] = spintbu ^ flip;

}

#define MAXWALKERS 6

__global__ void d_Walk(const int dir, int **neigAll[],
		       int n,
		       int nrep,
		       int nclone,
		       int nc,
		       const unsigned char *rotate,
		       MYWORD **walker0All[],
		       MYWORD **walker1All[]) {

  const unsigned int sampleId = blockIdx.y;

  const int **neig = (const int **) neigAll[sampleId];

  const unsigned int tid = threadIdx.x+blockDim.x*blockIdx.x;

  const unsigned int iclone = tid/(n*nrep);
  const unsigned int cloneioff = iclone*n*nrep;

  int tidModNNRep = tid - (iclone*n*nrep);

  int whichrep  = tidModNNRep % nrep;
  int whichsite = tidModNNRep / nrep;

  if(iclone>=nclone || whichsite>=n) return;

  MYWORD *walker0=dir ? walker0All[sampleId][1]:walker0All[sampleId][0];
  MYWORD *Owalker0=dir ? walker0All[sampleId][0]:walker0All[sampleId][1];  
  MYWORD *walker1=dir ? walker1All[sampleId][1]:walker1All[sampleId][0];
  MYWORD *Owalker1=dir ? walker1All[sampleId][0]:walker1All[sampleId][1];
  
  neig += DEGREE*dir;

  Owalker0 += cloneioff + whichrep;
  Owalker1 += cloneioff + whichrep;  
  unsigned char rot = rotate[whichsite];
  for(int ic=0; ic<nc; ic++) {
  
    MYWORD C0 = LDG(walker0[cloneioff + whichsite*nrep + whichrep]);
    MYWORD C1 = LDG(walker1[cloneioff + whichsite*nrep + whichrep]);  

    const int punto_n = neig[ic][whichsite];
    MYWORD scra0 = Owalker0[punto_n*nrep];
    MYWORD scra1 = Owalker1[punto_n*nrep];  

    switch(ic) {
    	case 0:
	   if (rot & 0x01) {
       	       scra0 = RotateE(scra0);
       	       scra1 = RotateE(scra1);
	       C0 = RotateO(C0);
       	       C1 = RotateO(C1);       
           }
	   break;
	 case 1:
	   if (rot & 0x02) {
	      scra0 = RotateN(scra0);
	      scra1 = RotateN(scra1);
	      C0 = RotateS(C0);
       	      C1 = RotateS(C1);       
	   }
	   break;
	 case 2:
	   if (rot & 0x04) {
	      scra0 = RotateU(scra0);
	      scra1 = RotateU(scra1);
	      C0 = RotateD(C0);
       	      C1 = RotateD(C1);       
	   }
	   break;
	 case 3:
	   if (rot & 0x08) {
	      scra0 = RotateO(scra0);
	      scra1 = RotateO(scra1);
	      C0 = RotateE(C0);
       	      C1 = RotateE(C1);       
	   }
	   break;
	 case 4:
	   if (rot & 0x10) {
	      scra0 = RotateS(scra0);
	      scra1 = RotateS(scra1);
	      C0 = RotateN(C0);
       	      C1 = RotateN(C1);        
	   }
	   break;
	 case 5:
           if (rot & 0x20) {
	      scra0 = RotateD(scra0);
	      scra1 = RotateD(scra1);
	      C0 = RotateU(C0);
       	      C1 = RotateU(C1);       
	   }
	   break;
	 default:
	   printf("Invalid walker's number %d, max is %d\n",ic,MAXWALKERS);
	   return;
    }
    walker0[cloneioff + whichsite*nrep + whichrep]=scra0;      
    walker1[cloneioff + whichsite*nrep + whichrep]=scra1;
    Owalker0[punto_n*nrep]=C0;
    Owalker1[punto_n*nrep]=C1;

    if(ic<(nc-1)){
      walker0 += HALF_MSC_VNRNBETAS;
      walker1 += HALF_MSC_VNRNBETAS;
      Owalker0 += HALF_MSC_VNRNBETAS;
      Owalker1 += HALF_MSC_VNRNBETAS;
    }
  }
}

#undef MASK_E
#undef MASK_N
#undef MASK_U
#undef MASK_O
#undef MASK_S
#undef MASK_D

#undef NMASK_E
#undef NMASK_N
#undef NMASK_U
#undef NMASK_O
#undef NMASK_S
#undef NMASK_D

// MMCUDA

#define MAKEMATR_RC 1
#if !defined(TRUE)
enum {FALSE, TRUE};
#endif
#if !defined(MAKEMATR_RC) 
#define MAKEMATR_RC 12
#endif

void **mmcuda(void ***rp, int r, int c, int s, int init) {
  int i;
  char **pc;
  short int **psi;
  int **pi;
  double **pd;
  char **d_pc;
  short int **d_psi;
  int **d_pi;
  double **d_pd;


  switch(s) {
  case sizeof(char):
    pc=(char **)malloc(r*sizeof(char *));
    if(!pc) create_error("error in makematr 1\n");
    MY_CUDA_CHECK( cudaMalloc( (void **) &d_pc, r*sizeof(char*) ) );
    for(i=0; i<r; i++) {
      MY_CUDA_CHECK( cudaMalloc( (void **) &pc[i], c*sizeof(char) ) );
      if(init) {
            MY_CUDA_CHECK( cudaMemset( pc[i], 0, c*sizeof(char) ) );
      }
    }
    MY_CUDA_CHECK( cudaMemcpy( d_pc, pc, r*sizeof(char *), cudaMemcpyHostToDevice ) );
    rp[0]=(void **)d_pc;
    return (void **)pc;
  case sizeof(short int):
    psi=(short int **)malloc(r*sizeof(short int*));
    if(!psi) create_error( "error in makematr 2\n");
    MY_CUDA_CHECK( cudaMalloc( (void **) &d_psi, r*sizeof(short int*) ) );
    for(i=0; i<r; i++) {
      MY_CUDA_CHECK( cudaMalloc( (void **) &psi[i], c*sizeof(short int) ) );
      if(init) {
            MY_CUDA_CHECK( cudaMemset( psi[i], 0, c*sizeof(short int) ) );
      }
    }
    MY_CUDA_CHECK( cudaMemcpy( d_psi, psi, r*sizeof(short int*), cudaMemcpyHostToDevice ) );
    rp[0]=(void **)d_psi;
    return (void **)psi;
  case sizeof(int):
    pi=(int **)malloc(r*sizeof(int*));
    if(!pi) create_error( "error in makematr 3\n");
    MY_CUDA_CHECK( cudaMalloc( (void **) &d_pi, r*sizeof(int*) ) );
    for(i=0; i<r; i++) {
      MY_CUDA_CHECK( cudaMalloc( (void **) &pi[i], c*sizeof(int) ) );
      if(init) {
            MY_CUDA_CHECK( cudaMemset( pi[i], 0, c*sizeof(int) ) );
      }
    }
    MY_CUDA_CHECK( cudaMemcpy( d_pi, pi, r*sizeof(int *), cudaMemcpyHostToDevice ) );
    rp[0]=(void **)d_pi;
    return (void **)pi;
  case sizeof(double):
    pd=(double **)malloc(r*sizeof(double*));
    if(!pd) create_error( "error in makematr 4 for %d rows\n",r);
    MY_CUDA_CHECK( cudaMalloc( (void **) &d_pd, r*sizeof(double*) ) );
    for(i=0; i<r; i++) {
      MY_CUDA_CHECK( cudaMalloc( (void **) &pd[i], c*sizeof(double) ) );
      if(init) {
            MY_CUDA_CHECK( cudaMemset( pd[i], 0, c*sizeof(double) ) );
      }
    }
    MY_CUDA_CHECK( cudaMemcpy( d_pd, pd, r*sizeof(double *), cudaMemcpyHostToDevice ) );
    rp[0]=(void **)d_pd;
    return (void **)pd;
  default:
    create_error("Unexpected size: %d\n",s);
    break;
  }
  return NULL;
}

#include "../include/header.h"
#include <signal.h>
#include <sys/stat.h>
#include <unistd.h>
#include <libgen.h>

char dir[MAXSTR];
char dirsample[MAXNSAMPLES][MAXSTR];
char dir_rep[MAXNSAMPLES][NR][MAXSTR];

char my_host_name[100];
int my_pid;

#define V8 (V>>3)
char janus_uu[NBETAS*V8];

void write_conf(unsigned long long ib, int nbits){
  char nameconf[MAXSTR], nameconf_t[MAXSTR];
  FILE *Fconfig;
  int ibeta, iclon, ir, ibit;

  struct stat statbuf;
  size_t predicted_size = sizeof(data) + sizeof(double)*NBETAS
    +sizeof(s_aleatorio_HQ_64bits)+sizeof(s_xoshiro256pp)
    +sizeof(uint8_t)*NBETAS + sizeof(char)*NBETAS*V8
    +sizeof(int)*NBETAS;

  int error=0;

  for(ibit=0;ibit<nbits;ibit++){
    unpacking_u(ibit);

    data.seed_J = seeds_J[ibit]; //seed of this sample
    
    for(ir=0;ir<NR;ir++){
      Janus_packing_for_write(ibit,ir);

      sprintf(nameconf_t,"%s/conf.%016llu",dir_rep[ibit][ir],ib);
      if( (Fconfig=fopen(nameconf_t,"wb"))==NULL ){
	create_error("Problems opening %s\n",nameconf_t);
	error=1;
	break;
      }
      
      fwrite(&data,sizeof(s_data), 1L, Fconfig);
      fwrite(betas,sizeof(double), NBETAS, Fconfig);

      fwrite(&random_PRC[ir], sizeof(s_aleatorio_HQ_64bits), 1L, Fconfig);
      fwrite(&random_xoshiro256pp[ir], sizeof(s_xoshiro256pp), 1L, Fconfig);

      fwrite(which_clon_this_beta[ibit][ir], sizeof(uint8_t), NBETAS, Fconfig);

      if( fwrite(janus_uu,sizeof(char), (size_t) NBETAS*V8, Fconfig) != (NBETAS*V8) ){
	create_error("Problems writing conf in %s\n",nameconf_t);
	error=1;
	break;
      }

      for(ibeta=0;ibeta<NBETAS;ibeta++){
	iclon = which_clon_this_beta[ibit][ir][ibeta];
	fwrite(&h_Ener[iclon+ir*NBETAS+ibit*NRNBETAS],sizeof(int), 1L, Fconfig);
      }
    
      fclose(Fconfig);

      sprintf(nameconf,"%s/conf",dir_rep[ibit][ir]);
      remove(nameconf);
      if( link(nameconf_t,nameconf) ){
	create_error("Problems in hard-lin to %s\n",nameconf_t);
	error=1;
	break;
      }

      stat(nameconf, &statbuf);
      if (statbuf.st_nlink!=2){
	create_error("The number of hard-link %s is not 2\n",nameconf);
	error=1;
	break;
      }

      if(statbuf.st_size != predicted_size ){
	create_error("There are only %d bytes in %s (must be %d)\n", statbuf.st_size, nameconf, predicted_size);
	error=1;
	break;
      }
    
    }

    if(error==1)
      break;
  }

  if(error==1)
    print_and_exit("Problems writing confs\n");
  
}

void read_conf(int nbits)
{
  char nameconf[MAXSTR];
  FILE *Fconfig;
  int ibeta, iclon, ir, ibit;
  s_data data_c;
  double betas_c[NBETAS];
  int Ener, Ener_scalar;
  
  struct stat statbuf;
  size_t predicted_size = sizeof(data) + sizeof(double)*NBETAS
    +sizeof(s_aleatorio_HQ_64bits)+sizeof(s_xoshiro256pp)
    +sizeof(uint8_t)*NBETAS + sizeof(char)*NBETAS*V8
    +sizeof(int)*NBETAS;

  int error=0;

  for(ibit=0;ibit<nbits;ibit++){
    for(ir=0;ir<NR;ir++){
      sprintf(nameconf,"%s/conf",dir_rep[ibit][ir]);

      stat(nameconf, &statbuf);
      if(statbuf.st_size != predicted_size ){
	create_error("There are only %d bytes in %s (must be %d)\n",
		     statbuf.st_size, nameconf, predicted_size);
	error=1;
	break;
      }

      if( (Fconfig=fopen(nameconf,"rb"))==NULL ){
	create_error("Problems opening %s\n",nameconf);
	error=1;
	break;
      }

      fread(&data_c,sizeof(s_data), 1L, Fconfig);
      check_data(data_c);
      if ( data.nc != data_c.nc ||
	   data.nterm != data_c.nterm){
	create_error("Different data from input and %s\n",nameconf);
	error=1;
	break;
      }

      if(seeds_J[ibit] != data_c.seed_J){
	create_error("Different seed_J generated and read in %s\n", nameconf);
	error=1;
	break;
      }
      
      data.now = data_c.now;
      stop(&data.maxtime);

      fread(betas_c,sizeof(double), NBETAS, Fconfig);
      for(ibeta=0;ibeta<NBETAS;ibeta++)
	if(fabs(betas_c[ibeta]-betas[ibeta])> 1e-12){
	  create_error("Different beta (ibeta=%d, %.12g != %.12g)\n",
		       ibeta, betas_c[ibeta],betas[ibeta]);
	  error=1;
	  break;
	}

      fread(&random_PRC[ir], sizeof(s_aleatorio_HQ_64bits), 1L, Fconfig);
      fread(&random_xoshiro256pp[ir], sizeof(s_xoshiro256pp), 1L, Fconfig);
      
      fread(which_clon_this_beta[ibit][ir], sizeof(uint8_t), NBETAS, Fconfig);

      for(iclon=0;iclon<NBETAS;iclon++)
	which_beta_this_clon[ibit][ir][iclon]=NBETAS;

      for(ibeta=0;ibeta<NBETAS;ibeta++)
	which_beta_this_clon[ibit][ir][which_clon_this_beta[ibit][ir][ibeta]] = ibeta;

      for(iclon=0;iclon<NBETAS;iclon++)
	if(which_beta_this_clon[ibit][ir][iclon] == NBETAS ){
	  create_error("There is an error in the permutatuons in %s\n",nameconf);
	  error = 1;
	  break;
	}
      
      if( fread(janus_uu,sizeof(char), (size_t) NBETAS*V8, Fconfig) != (NBETAS*V8) ){
	create_error("Problems reading conf in %s\n",nameconf);
	error=1;
	break;
      }

      Janus_unpacking_for_read(ibit,ir);
      
      for(ibeta=0;ibeta<NBETAS;ibeta++){
	fread(&Ener,sizeof(int), 1L, Fconfig);
	iclon = which_clon_this_beta[ibit][ir][ibeta];
	Ener_scalar=calculate_scalar_energy(ibit,ir,iclon);
	if ( Ener != Ener_scalar ){
	  create_error("Energy read and calculate differ in ibit=%d ir=%d ibeta=%d:escalar=%d, read=%d\n",
		       ibit, ir, ibeta, Ener_scalar, Ener);
	  error=1;
	  break;
	}
      }

      fclose(Fconfig);
    
      if(error==1)
	break;
    }

    if(error==1)
      print_and_exit("Problems reading confs\n");
  }

  packing_u(nbits);
  check_packing(nbits);
}

void Janus_packing_for_write(int ibit, int ir)
{
  int site, bit, x, y, z, iclon;
  int pos;
  
  memset(janus_uu,0,sizeof(char)*NBETAS*V8);

  for(iclon=0;iclon<NBETAS;iclon++){
    site = 0;
    for(z=0;z<Lz;z++)
      for(y=0;y<L;y++)
	for(x=0;x<L;x++){
	  bit = site&7;
	  pos = (site>>3) + iclon*V8;
	  janus_uu[pos] |= uu[ibit][ir][iclon][site]<<bit;
	  site++;
	}
  }
  
}

void Janus_unpacking_for_read(int ibit, int ir)
{
  int site, bit, iclon, site_MSC;
  int pos;
  
  pos = 0;
  for(iclon=0;iclon<NBETAS;iclon++){
    site = 0;
    for(site_MSC=0;site_MSC<V8;site_MSC++){
      for(bit=0;bit<8;bit++){
	uu[ibit][ir][iclon][site] = (janus_uu[pos]>>bit)&1;
	site++;
      }
      pos++;
    }
  }

}

#undef V8

void write_walkers(unsigned long long ib, int nbits, int prev_ib){
  char nameconf[MAXSTR], nameconf_t[MAXSTR];
  FILE *Fconfig;
  int ibit;

  struct stat statbuf;
  size_t predicted_size = sizeof(icds) +
    +NR*sizeof(s_aleatorio_HQ_64bits)+NR*sizeof(s_xoshiro256pp)
    +4*sizeof(MYWORD)*HALF_MSC_VNRNBETAS*data.nc;

  int error = 0;
  for(ibit=0;ibit<nbits;ibit++){
    sprintf(nameconf_t,"%s/walkers.%016llu",dirsample[ibit],ib);
    if( (Fconfig=fopen(nameconf_t,"wb"))==NULL ){
      create_error("Problems opening %s\n",nameconf_t);
      error=1;
      break;
    }
    fwrite(p2icds[ibit],sizeof(icds), 1L, Fconfig);
    fwrite(&write_random_PRC_C[ibit], sizeof(s_aleatorio_HQ_64bits), (size_t) NR, Fconfig);
    fwrite(&write_random_xoshiro256pp_C[ibit],sizeof(s_xoshiro256pp), (size_t) NR, Fconfig);

    if( fwrite(write_CPUwalker0[ibit], sizeof(MYWORD),
	       (size_t) 2*HALF_MSC_VNRNBETAS*data.nc, Fconfig) != (2*HALF_MSC_VNRNBETAS*data.nc) ){
      create_error("Problems writing walkers in %s\n",nameconf_t);
      error=1;
      break;
    }

    if( fwrite(write_CPUwalker1[ibit], sizeof(MYWORD),
	       (size_t) 2*HALF_MSC_VNRNBETAS*data.nc, Fconfig) != (2*HALF_MSC_VNRNBETAS*data.nc) ){
      create_error("Problems writing walkers in %s\n",nameconf_t);
      error=1;
      break;
    }
    
    fclose(Fconfig);
    
    sprintf(nameconf,"%s/walkers",dirsample[ibit]);
    remove(nameconf);
    if( link(nameconf_t,nameconf) ){
      create_error("Problems in hard-lin to %s\n",nameconf_t);
      error=1;
      break;
    }

    stat(nameconf, &statbuf);
    if (statbuf.st_nlink!=2){
      create_error("The number of hard-link %s is not 2\n",nameconf);
      error=1;
      break;
    }

    if(statbuf.st_size != predicted_size ){
      create_error("There are only %d bytes in %s (must be %d)\n", statbuf.st_size, nameconf, predicted_size);
      error=1;
      break;
    }

    if(prev_ib){
      sprintf(nameconf_t,"%s/walkers.%016llu",dirsample[ibit],prev_ib);
      remove(nameconf_t);
    }
    
  }
  
  if(error==1)
    print_and_exit("Problems writing walkers\n");
  
}

void read_walkers(int nbits)
{

  char nameconf[MAXSTR];
  FILE *Fconfig;
  int ibit;
  icds p2icds_c;
  int error = 0;
  
  struct stat statbuf;
  size_t predicted_size = sizeof(icds) +
    +NR*sizeof(s_aleatorio_HQ_64bits)+NR*sizeof(s_xoshiro256pp)
    +4*sizeof(MYWORD)*HALF_MSC_VNRNBETAS*data.nc;

  for(ibit=0;ibit<nbits;ibit++){
    sprintf(nameconf,"%s/walkers",dirsample[ibit]);
    
    stat(nameconf, &statbuf);
    if(statbuf.st_size != predicted_size ){
      create_error("There are only %d bytes in %s (must be %d)\n",
		   statbuf.st_size, nameconf, predicted_size);
      error=1;
      break;
    }
  
    if( (Fconfig=fopen(nameconf,"rb"))==NULL ){
      create_error("Problems opening %s\n",nameconf);
      error=1;
      break;
    }

    fread(&p2icds_c,sizeof(icds), 1L, Fconfig);
    if( (p2icds[ibit]->isample != p2icds_c.isample) ||
	(p2icds[ibit]->Nrep != p2icds_c.Nrep) ||
	(p2icds[ibit]->NUM_WALKERS_PER_SITE != p2icds_c.NUM_WALKERS_PER_SITE) ||
	(p2icds[ibit]->Lv != p2icds_c.Lv) ||
	(p2icds[ibit]->seed != p2icds_c.seed) ||
	(fabs(p2icds[ibit]->beta - p2icds_c.beta)>1e-6) ){
      create_error("Different icds from simulation and %d\n", nameconf);
      error = 1;
      break;
    }
    
    p2icds[ibit]->first = p2icds_c.first;
 
    fread(&random_PRC_C[ibit], sizeof(s_aleatorio_HQ_64bits), (size_t) NR, Fconfig);
    fread(&random_xoshiro256pp_C[ibit],sizeof(s_xoshiro256pp), (size_t) NR, Fconfig);

    if( fread(CPUwalker0[ibit], sizeof(MYWORD),
	      (size_t) 2*HALF_MSC_VNRNBETAS*data.nc, Fconfig) != (2*HALF_MSC_VNRNBETAS*data.nc) ){
      create_error("Problems reading walkerss in %s\n",nameconf);
      error=1;
      break;
    }
  
    if( fread(CPUwalker1[ibit], sizeof(MYWORD),
	      (size_t) 2*HALF_MSC_VNRNBETAS*data.nc, Fconfig) != (2*HALF_MSC_VNRNBETAS*data.nc) ){
      create_error("Problems writing walkers in %s\n",nameconf);
      error=1;
      break;
    }
    
    fclose(Fconfig);
  }

  if(error==1)
    print_and_exit("Problems reading walkers\n");
  
}

void write_measures(unsigned long long ib, int nbits)
{
  FILE *Fmeasures;
  char namefile[MAXSTR];

  int ibeta, ibit;

  double Bc;
  
  //PRINT in STDOUT
  printf("Block = %d\n", ib);
  for(ibit=0;ibit<nbits;ibit++){
    printf("Sample = %d\n",list_samples[ibit]);
    for(ibeta=0;ibeta<NBETAS;ibeta++){
      Bc = 0.25*log(1.+1./average_Energy_C[ibit][ibeta]);
      printf("ibeta=%d: E_s=%10.8g E_c=%10.8g T_cl=%10.8g ",
	     ibeta, average_Energy[ibit][ibeta], average_Energy_C[ibit][ibeta],1./Bc);
    }
    printf("\n");
  }
  fflush(stdout);

  //MEASURES
  for(ibit=0;ibit<nbits;ibit++){
    sprintf(namefile,"%sevol_E_spins.txt",dirsample[ibit]);
    if ( (Fmeasures=fopen(namefile,"a"))==NULL ){
      create_error("Problems opening %s\n", namefile);
      print_and_exit("Problems writing measures\n");
    }
    
    fprintf(Fmeasures, "%llu ", ib);
    for(ibeta=0;ibeta<NBETAS;ibeta++)
      fprintf(Fmeasures, "%14.10g %14.10g", average_Energy[ibit][ibeta], std_Energy[ibit][ibeta]);
    
    fprintf(Fmeasures,"\n");
    fclose(Fmeasures);

    sprintf(namefile,"%sevol_E_walkers.txt",dirsample[ibit]);
    if ( (Fmeasures=fopen(namefile,"a"))==NULL ){
      create_error("Problems opening %s\n", namefile);
      print_and_exit("Problems writing measures\n");
    }
    
    fprintf(Fmeasures, "%llu ", ib);
    for(ibeta=0;ibeta<NBETAS;ibeta++)
      fprintf(Fmeasures, "%14.10g %14.10g", average_Energy_C[ibit][ibeta], std_Energy_C[ibit][ibeta]);
    
    fprintf(Fmeasures,"\n");
    fclose(Fmeasures);

    sprintf(namefile,"%sevol_T_walkers.txt",dirsample[ibit]);
    if ( (Fmeasures=fopen(namefile,"a"))==NULL ){
      create_error("Problems opening %s\n", namefile);
      print_and_exit("Problems writing measures\n");
    }
    
    fprintf(Fmeasures, "%llu ", ib);
    for(ibeta=0;ibeta<NBETAS;ibeta++){
      Bc = 0.25*log(1.+1./average_Energy_C[ibit][ibeta]);
      fprintf(Fmeasures, "%14.10g ", 1./Bc);
    }
    fprintf(Fmeasures,"\n");
    fclose(Fmeasures);
    
    
  }//ibit
  
}

void read_input(const char *name){
  FILE *Finput;
  int j;
  unsigned long long *ptdata_ull;
  int *ptdata_int;
  randint *ptdata_randint;
  char line[MAXSTR+1];

  if( NULL==(Finput=fopen(name,"r")) )
    print_and_exit("Error: No file %s\n",name);
  
  fgets(line,MAXSTR,Finput);
  sscanf(line,"%s",dir);

  for (j=0,ptdata_ull=&data.maxtime;j<NDATULL;j++){
    fgets(line,MAXSTR,Finput);
    sscanf(line,"%llu",ptdata_ull++);
  }
  
  for (j=0,ptdata_int=&data.ntimes;j<NDATINT;j++){
    fgets(line,MAXSTR,Finput);
    sscanf(line,"%u",ptdata_int++);
  }

  for (j=0,ptdata_randint=&data.seed_J;j<NDATRAND;j++){
    fgets(line,MAXSTR,Finput);
    sscanf(line,"%llu",ptdata_randint++);
  }

  fclose(Finput);

  check_data(data);

  printf("Path: %s\n",dir);

  print_data(stdout, &data);
  fflush(stdout);

}

void check_data(s_data d)
{
  if (d.maxtime<0 || d.now<0 || d.nterm<0 || d.nc < 0|| d.ntimes<0)
    print_and_exit("Bad data input\n");
    
  // l= M<<16 + L
  if ( d.nr != NR)
    print_and_exit("d.nr = %d != NR = %d\n", d.nr, NR);
      
  if (d.l!=L)
    print_and_exit("l=%d != L=%d\n",d.l,L);

}

void print_data(FILE *file, s_data *dat){
    fprintf(file,"maxtime      %16llu\n",dat->maxtime);
    fprintf(file,"now          %16llu\n",dat->now);
    fprintf(file,"ntimes       %16d\n",dat->ntimes);
    fprintf(file,"nterm        %16d\n",dat->nterm);
    fprintf(file,"nr           %16d\n",dat->nr);
    fprintf(file,"nc           %16d\n",dat->nc);
    fprintf(file,"L            %16d\n",((dat->l)&0xffff));
    fprintf(file,"M            %16d\n",((dat->l)>>16));
    fprintf(file,"flag         %16d\n",dat->flag);
    fprintf(file,"seed_J       %016llx\n",dat->seed_J);
    fprintf(file,"seed_u       %016llx\n",dat->seed_u);
    for(int ir=0;ir<NR;ir++)
      fprintf(file,"seed_MC[%d]   %016llx\n",ir,dat->seed_MC[ir]);
}

void read_betas(const char * name_betas)
{
  int ibeta;

  double dummy;
  FILE *fi_betas;
  
  if (NULL==(fi_betas=fopen(name_betas,"rt")))
    print_and_exit("File %s doesn't exist.\n",name_betas);

  for (ibeta=0;ibeta<NBETAS;ibeta++)
    if (EOF==fscanf(fi_betas,"%lf ",&betas[ibeta]))
      print_and_exit("File %s only has %d lines (it needs %d)\n",
		     name_betas,ibeta,NBETAS);
  // More number of betas in file?
  if (EOF!=fscanf(fi_betas,"%lf",&dummy))
    fprintf(stderr,"More lines in file %s\n",name_betas); 
  
  for(ibeta=0;ibeta<NBETAS-1;ibeta++){
    if (betas[ibeta]<betas[ibeta+1])
      print_and_exit("Error in file %s: betas don't sort from largest to smallest!\n",name_betas);
  }

  fclose(fi_betas);
}


void read_lut(const char * name, const char * name_betas)  
{
  FILE *Fin;
  double betas_read[NBETAS];
  int ib;

  if(NULL==(Fin=fopen(name,"rb")))
    print_and_exit("Problems opening %s\n",name);

  for(ib=0;ib<NBETAS;ib++)
    betas_read[ib]=0;
  
  if(1!=fread(betas_read,sizeof(betas_read),(size_t) 1,Fin)){
    printf("betas_read[0]=%.14g\n",betas_read[0]);
    printf("betas_read[1]=%.14g\n",betas_read[1]);
    print_and_exit("Problems reading betas\n");
  }

  for(ib=0;ib<NBETAS;ib++){
    if(fabs(betas[ib]-betas_read[ib])>1e-14)
      print_and_exit("betas[%d]=%.14g from %s, betas[%d]=%.14g from %s (difference=%.14g)\n",
                     ib,betas[ib],name_betas,ib,betas_read[ib],name,betas[ib]-betas_read[ib]);
  }


  if(1!=fread(h_LUT,sizeof(s_lut_heat_bath)*NBETAS,(size_t) 1,Fin))
    print_and_exit("Problems reading lut_heat_bath\n");
  
  fclose(Fin);
  
}

void read_list_samples(const char * name, int nbits)
{
  int ibit=0;
  int dummy;
  FILE *Fin;
  
  if (name[0] == '\0'){
    for(ibit=0;ibit<nbits;ibit++)
      list_samples[ibit] = ibit;
  }else{
  
    if (NULL==(Fin=fopen(name,"rt")))
      print_and_exit("File %s doesn't exist.\n",name);

    for (ibit=0;ibit<nbits;ibit++)
      if (EOF==fscanf(Fin,"%d ",&list_samples[ibit]))
	print_and_exit("File %s only has %d lines (it needs %d)\n",
		       name,ibit,nbits);
    // More number of betas in file?
    if (EOF!=fscanf(Fin,"%d",&dummy))
      fprintf(stderr,"More lines in file %s\n",name); 
    
    fclose(Fin); 
  }

  for(ibit=0;ibit<(nbits-1);ibit++)
    if(list_samples[ibit]>=list_samples[ibit+1]){
      create_error("List of samples %s is not sorted.\n", name);
      print_and_exit("List of samples must be sorted!\n");
    }
  
}

void read_list_times(const char * name)
{
  int itime;
  unsigned long long i;
  unsigned long long dummy;
  FILE *Fin;

  if(NULL==(list_times = (unsigned long long *)malloc(sizeof(unsigned long long)*data.ntimes)))
    print_and_exit("Problems allocating memory for list_times\n");
  
  if (NULL==(Fin=fopen(name,"rt")))
    print_and_exit("File %s doesn't exist.\n",name);
  
  for(itime=0;itime<data.ntimes;itime++)
    if (EOF==fscanf(Fin,"%llu ",&list_times[itime]))
      print_and_exit("File %s only has %d lines (it needs %d)\n",
		     name,itime,data.ntimes);
  
  if ( EOF!=fscanf(Fin,"%llu",&dummy) ){
    fprintf(stderr,"More lines in file %s\n",name); 
  }
  
  fclose(Fin);
  
  for(itime=0;itime<(data.ntimes-1);itime++)
    if(list_times[itime]>=list_times[itime+1]){
      create_error("List of times %s is not sorted.\n", name);
      print_and_exit("List of times must be sorted!\n");
    }
}

void create_sample_path(int sample, int nbits)
{
  int ir, ibit;
  
  sprintf(dir,"%sI%03d/",dir,sample);
  mkdir(dir,(mode_t) 0775);

  for(ibit=0;ibit<nbits;ibit++){
    sprintf(dirsample[ibit],"%sBIT%03d/",dir,list_samples[ibit]);
    mkdir(dirsample[ibit],(mode_t) 0775);

    for(ir=0;ir<NR;ir++){
      sprintf(dir_rep[ibit][ir],"%sR%02d/",dirsample[ibit],ir);
      mkdir(dir_rep[ibit][ir], (mode_t) 0775);
    }
  }
  
  printf("Dirsample = %s\n",dir);
  fflush(stdout);

}

/* The possible values of flag for the simulation are:

   flag = 0 ----> Backup
   flag = 1 ---> new run
   flag = 2 ---> Backup if there are conf, new run in other case

   flag = -1 ---> New run from CUBE packing program
   flag = -2 ---> Backup reading couplings from CUBE packing programa
   
*/
void check_and_prepare_simulation_backup(int sample, int nbits)
{

  FILE *Finput;
  int there_is_conf, check_now = 0;
  int ir, ibit, now = 0;
  char nameconf[MAXSTR];
  s_data data_c;

  for(ibit=0; ibit<nbits;ibit++){
    there_is_conf=0;
    for(ir=0;ir<NR;ir++){
      sprintf(nameconf,"%s/conf",dir_rep[ibit][ir]);
      if ((Finput=fopen(nameconf,"r"))!=NULL){
	fread(&data_c,sizeof(s_data), 1L, Finput);
	if((check_now>0) && (data_c.now != now)){
	  create_error("Different number of iterations in samples from %s\n", dir);
	  print_and_exit("Problems with conf files\n");
	}
	now = data_c.now;
	there_is_conf++;
	check_now++;
	fclose(Finput);
      }
    }

    //Check how many conf files are there?
    if (there_is_conf>0){
      if(there_is_conf!=NR){
	create_error("There are less conf files than it must be in %s (%d/%d)\n",
		     dirsample[ibit],there_is_conf,NR);
	print_and_exit("Problems with conf files\n");
      }
    }
  }
  
  if(data.flag==2){
    if(there_is_conf){
      data.flag = 0;
    }else{
      data.flag = 1;
    }
  }
  
  if ( data.flag == 1 ){

    if(there_is_conf){
      create_error("There are some confs in %s.\nIf you want to make a new run, please delete the configuration manually\n", dir);
      print_and_exit("There are confs\n");
    }
    
  }else if ( data.flag == 0 ){
    if(there_is_conf==0){
      create_error("There is not %s, it is impossible to do a backup\n",
		     nameconf);
      print_and_exit("Impossible to do a backup\n");
    }
  }else{
    print_and_exit("data.flag is not recognize\n");
  }

  write_seeds = get_seeds(nbits);
  
}

int get_seeds(int nbits)
{
  
  char namefile[MAXSTR];
  FILE *c_file;
  unsigned long long sJ, su, sMC[NR];
  int n_fields;
  int no_seeds = 0;
  int ir, ibit;
  uint64_t temp_seeds[MAXNSAMPLES];
  
  sprintf(namefile,"%s/master_seedJ.in", dir);

  if ( (c_file=fopen(namefile,"r"))==NULL ){
    if (data.flag == 1){
      no_seeds++;
    }
    else{
      sleep(1);
      if ( (c_file=fopen(namefile,"r"))==NULL ){
	create_error("There is not file %s, and it must exist\n", namefile);
	print_and_exit("Problems reading seeds\n");
	}
    }
  }

  if(no_seeds==0){
    if(EOF==fscanf(c_file,"%llu",&sJ))
      print_and_exit("Nothing in file %s\n",namefile);
    
    fclose(c_file);
    
    if( (data.seed_J>0) && (data.seed_J != sJ) ){
      create_error("Master seedJ in input file and %s differ\n",namefile);
      print_and_exit("Problems reading master seedJ\n");
    }else{
      data.seed_J = sJ;
    }
  }

  generate_seeds_from_one(&data.seed_J, temp_seeds, MAXNSAMPLES);
  for(ibit=0; ibit<nbits; ibit++){
    seeds_J[ibit] = (randint) temp_seeds[list_samples[ibit]];
#ifdef JANUS
    seeds_J[ibit] = 9164166177945067407LLU;
#endif
  }

  memset(countJ0,0, sizeof(int)*MAXNSAMPLES);
  
  for(ibit=0; ibit<nbits;ibit++){

    n_fields = 0;
    sprintf(namefile,"%s/seeds.in", dirsample[ibit]);

    if ( (c_file=fopen(namefile,"r"))==NULL ){
      if (data.flag == 1){
	 no_seeds++;
	 continue;
      }
      else{
	sleep(1);
	if ( (c_file=fopen(namefile,"r"))==NULL ){
	  create_error("There is not file %s, and it must exist\n", namefile);
	  print_and_exit("Problems reading seeds\n");
	}
      }
    }
  
    n_fields += fscanf(c_file,"%llu %llu",&sJ,&su);
    for(ir=0;ir<NR;ir++)
      n_fields += fscanf(c_file," %llu", &sMC[ir]);

    n_fields += fscanf(c_file," %d", &countJ0[ibit]);

    if (n_fields!=(NR+3)){
      create_error("File %s must have 7 fields\n", namefile);
      print_and_exit("Problems reading seeds from %s\n",namefile);
      
    }
  
    fclose(c_file);

    //CHECKs

    if (seeds_J[ibit] != sJ) {
      create_error("Seed J generated and %s differ (bit=%d)\n", namefile, list_samples[ibit]);
      print_and_exit("Problems reading seeds\n");
    }else{
      seeds_J[ibit] = sJ;
    }
    
    if ( (data.seed_u>0) && (data.seed_u != su) ){
      create_error("Seed u in input and %s differ\n", namefile);
      print_and_exit("Problems reading seeds\n");
    }else{
      data.seed_u = su;
    }

    for(ir=0;ir<NR;ir++){
      if ( (data.seed_MC[ir]>0) && (data.seed_MC[ir] != sMC[ir]) ){
	create_error("Seed MC in input and %s differ\n", namefile);
	print_and_exit("Problems reading seeds\n");
      }else{
	data.seed_MC[ir] = sMC[ir];
      }
    }
  }
  
  if(no_seeds==(nbits+1)){ //+1 because the master_seedJ.in file
    return 1;
  }else if (no_seeds>0){
    print_and_exit("Some samples have file seeds.in, others not have\n");
  }
	   
  return 0;
}

void write_seeds_in_file(int nbits)
{
  char namefile[MAXSTR];
  FILE *Fseed;
  int ir, ibit;

  //WRITE MASTER SEED for J:
  sprintf(namefile,"%s/master_seedJ.in", dir);
  
  ignore_sigterm();
  if ( (Fseed=fopen(namefile,"w"))==NULL ){
    create_error("Impossible to open seed file to write\n");
    print_and_exit("Problems writing seeds\n");
  }

  fprintf(Fseed,"%llu\n", data.seed_J);
  fclose(Fseed);
    
  handle_sigterm();
  
  for(ibit=0;ibit<nbits;ibit++){
    sprintf(namefile,"%s/seeds.in", dirsample[ibit]);

    ignore_sigterm();
    if ( (Fseed=fopen(namefile,"w"))==NULL ){
      create_error("Impossible to open seed file to write\n");
      print_and_exit("Problems writing seeds\n");
    }

    fprintf(Fseed,"%llu %llu", seeds_J[ibit], data.seed_u);
    for( ir=0; ir<NR; ir++)
      fprintf(Fseed," %llu", data.seed_MC[ir]);

    fprintf(Fseed," %d\n",countJ0[ibit]);
    
    fclose(Fseed);
    
    handle_sigterm();
  }
}

void backup(int nbits)
{
  char error[MAXSTR];
  read_conf(nbits);
  read_walkers(nbits);

  print_data(stdout,&data);
  fflush(stdout);

  if (data.now >= data.maxtime){
    delete_running();
    create_error("There are %d iterations in %s\n", data.now, dirsample);
    print_and_exit("Problems in backup\n");
  }
  
}

void stop(unsigned long long int *newntimes)
{
  char name[MAXSTR];
  FILE *Fstop;
  
  unsigned long long *temp, dummy;
  int i, itime;
  
  sprintf(name,"%schangentimes.in",dir);
  Fstop=fopen(name,"rb");
  if (Fstop!=NULL){
    if (fscanf(Fstop,"%llu",temp)==EOF){
      fclose(Fstop);
      remove(name); // without this the run is blocked
      delete_running();
      printf(" File '%s' exists .\a\n",name);
    }
    *newntimes = (*temp);

    free(list_times);
    if(NULL==(list_times = (unsigned long long *)malloc(sizeof(unsigned long long)*(*newntimes))))
      print_and_exit("Problems allocating memory for list_times\n");

    for(itime=0;itime<(*newntimes);itime++){
      if (EOF==fscanf(Fstop,"%llu ",&list_times[itime]))
	print_and_exit("File %s only has %d lines (it needs %d)\n",
		       name,itime,(*newntimes));
      
      for(itime=0;itime<(data.ntimes-1);itime++)
	if(list_times[itime]>=list_times[itime+1]){
	  create_error("List of times %s is not sorted.\n", name);
	  print_and_exit("List of times must be sorted!\n");
	}
    }

    fclose(Fstop);
  }
  
}

void ignore_sigterm(void)
{
    signal(SIGTERM,SIG_IGN);
    signal(SIGINT,SIG_IGN);
    signal(SIGUSR1,SIG_IGN);
}

void handle_sigterm(void)
{
    signal(SIGTERM,stop_run);
    signal(SIGINT,stop_run);
    signal(SIGUSR1,stop_run);
}

void stop_run(int signum)
{
  if (signum==SIGTERM)    
    create_error("Received SIGTERM\n");  
  if (signum==SIGINT)    
    create_error("Received SIGINT\n");  
  if (signum==SIGUSR1)    
    create_error("Received SIGUSR1\n");  

  print_and_exit("Aborting programma (SIGNAL=%d)\n",signum);  
}

void create_error(const char *format, ...)
{
  va_list list;
  FILE *Foutput;
  char name[MAXSTR];

  snprintf(name,(size_t) MAXSTR,"%s/error",dir);
  Foutput=fopen(name,"w");

  va_start(list,format);
  vfprintf(Foutput,format,list);
  va_end(list);
  
  print_data(Foutput,&data);          

  print_data(stderr,&data);     

  delete_running();
}

void create_running(void) // si existe running se para
{
    FILE *Frunning,*Flog;
    char nombre[MAXSTR];
    char nombre_link[MAXSTR];
   
    gethostname(my_host_name,99L);// hostname es global
    my_pid=getpid();          // pid es global

    sprintf(nombre,"%srunning.%s.%d",dir,my_host_name,my_pid);  
    sprintf(nombre_link,"%srunning",dir);  

    if (NULL==(Frunning=fopen(nombre,"w")))
	print_and_exit("Unable to open %s\n",nombre);
    writelog(Frunning,time(NULL));
    fclose(Frunning);

    // Ahora hace un link a 'running'
    if (link(nombre,nombre_link)==-1){// Existe 'running'
	remove(nombre); // comentar para dejar rastro de los intentos
	print_and_exit("%s exists\n",nombre_link);
    }
// ya ningun otro proceso puede trabajar con esta sample
    remove(nombre);

    sprintf(nombre,"%srun.log",dir);
    if (NULL==(Flog=fopen(nombre,"a")))
	print_and_exit("Unable to open %s\n",nombre);
    writelog(Flog,time(NULL));
    fclose(Flog);

}


void renew_running(void) // running debe existir, con iguales 
{                            // x, y, board, hostname pid
    FILE *Frunning;
    char nombre[MAXSTR];
    char newhostname[100];
    int newpid;

    sprintf(nombre,"%s/running",dir);  

    if (NULL==(Frunning=fopen(nombre,"r+"))){
	create_error("File %s doesn't exist\n",nombre);
	print_and_exit("File %s doesn't exist\n",nombre);
    }
    
    fscanf(Frunning,"%s %d",newhostname,&newpid);	
    if (strcmp(my_host_name,newhostname)){
	create_error("Host has changed from %s to %s\n",
		       my_host_name,newhostname);
	print_and_exit("Host has changed from %s to %s\n",
		       my_host_name,newhostname);
    }
    
    if (my_pid!=newpid){
      create_error("PID has changed from %d to %d\n",my_pid,newpid);
      print_and_exit("PID has changed from %d to %d\n",my_pid,newpid);
    }
    
    rewind(Frunning); // se prepara para sobreescribir
    writelog(Frunning,time(NULL));
    fclose(Frunning);
}

void delete_running(void)
{
    char nombre[MAXSTR];
   
    sprintf(nombre,"%s/running",dir);
    remove(nombre);
}

void writelog(FILE *Fout, time_t t)
{
    fprintf(Fout,"%s %5d %10ld %s",my_host_name,my_pid,t,ctime(&t));	
}

void print_help(const char *progname)
{
  printf("Usage: %s isample ibit beta.dat input.in LUT list_times", progname);
  printf(" device [max_time [list_samples]]\n");
  printf("\nArguments:\n");
  printf("  isample              Sample index (int)\n");
  printf("  nbits                Number of subsample (int)\n");
  printf("  beta.dat             Path to beta.dat file (string)\n");
  printf("  input.in             Path to input.in file (string)\n");
  printf("  LUT                  Path to LUT file (string)\n");
  printf("  list_times           Path to the list_of_times file (string)\n");
  printf("  device               Device number (int)\n");
  printf("  max_time (opt.)      Maximum time allowed (unsigned long long)\n");
  printf("                       A 0 value disable this option.\n");
  printf("  list_samples (opt.)  Path to the list_of_samples file (string)\n");
  printf("\nUse -h to display this help message.\n");
  exit(EXIT_SUCCESS);
}

void print_and_exit(const char *format, ...)
{
  va_list list;
    
  va_start(list,format);
  vfprintf(stderr,format,list);
  va_end(list);
  exit(1);
}

#ifdef JANUS
int compact_dump(char* outdata,FILE* ostr,int num) {

    int is,ib,j,nwrote;
    char* buffer;

    buffer=(char *) calloc (num/8,sizeof(char));

    for (is=0,ib=0;is<num;is+=8,ib++) {
        buffer[ib]=0;
        for(j=0;j<8;j++) {
            buffer[ib]+=(outdata[is+j])<<j ;
        }
    }

    nwrote=fwrite(buffer,sizeof(char),num/8,ostr);

    free(buffer);

    return 8*nwrote;

}

void write_conf_Js(int nbits)
{
    FILE *fhs;
    char nombre_hs[1024];

    for(int ibit=0;ibit<nbits;ibit++){
      sprintf(nombre_hs, "%s/conf_Js", dirsample[ibit]);
      if((fhs=fopen(nombre_hs,"w"))==NULL)
          print_and_exit("No pude abrir el fichero de configs. de Js.\n");
      
      compact_dump(Jx[ibit],fhs,V);
      compact_dump(Jy[ibit],fhs,V);
      compact_dump(Jz[ibit],fhs,V);

      fclose(fhs);
    }

} 
#endif

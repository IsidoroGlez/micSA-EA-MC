#include "bits-random.h"

//IO

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

  fclose(fi_betas);
}


void lee_lut(char * name_betas)  
{
  FILE *Fin;
  char name[1024];
  double betas_read[NUMBETAS];
  int ib;
  
  sprintf(name,"LUT_for_PRNG_nbits%02d_NB%d.bin",NUMBITSPREBUSQUEDAS,NUMBETAS);
  if(NULL==(Fin=fopen(name,"rb")))
    print_and_exit("Problemas abriendo %s\n",name);

  for(ib=0;ib<NUMBETAS;ib++)
    betas_read[ib]=0;
  
  if(1!=fread(betas_read,sizeof(betas_read),(size_t) 1,Fin)){
    printf("betas_read[0]=%.14g\n",betas_read[0]);
    printf("betas_read[1]=%.14g\n",betas_read[1]);
    print_and_exit("problemas leyendo betas\n");
  }
  if(1!=fread(lut_heat_bath,sizeof(lut_heat_bath),(size_t) 1,Fin))
    print_and_exit("problemas leyendo lut_heat_bath\n");

  for(ib=0;ib<NUMBETAS;ib++){
    if(fabs(betas[ib]-betas_read[ib])>1e-14)
      print_and_exit("betas[%d]=%.14g from %s, betas[%d]=%.14g from %s (difference=%.14g)\n",
		     ib,betas[ib],name_betas,ib,betas_read[ib],name,betas[ib]-betas_read[ib]);
  }
  
  fclose(Fin);
}


void print_and_exit(const char *format, ...)
{
  va_list list;
    
  va_start(list,format);
  vfprintf(stderr,format,list);
  va_end(list);
  exit(1);
}


//Rutina tomada de "GNU radio"
void obten_crc(uint32_t *crc_thread_sums)
{

  static const uint32_t table_crc[256] = {
    0x00000000U,0x04C11DB7U,0x09823B6EU,0x0D4326D9U,
    0x130476DCU,0x17C56B6BU,0x1A864DB2U,0x1E475005U,
    0x2608EDB8U,0x22C9F00FU,0x2F8AD6D6U,0x2B4BCB61U,
    0x350C9B64U,0x31CD86D3U,0x3C8EA00AU,0x384FBDBDU,
    0x4C11DB70U,0x48D0C6C7U,0x4593E01EU,0x4152FDA9U,
    0x5F15ADACU,0x5BD4B01BU,0x569796C2U,0x52568B75U,
    0x6A1936C8U,0x6ED82B7FU,0x639B0DA6U,0x675A1011U,
    0x791D4014U,0x7DDC5DA3U,0x709F7B7AU,0x745E66CDU,
    0x9823B6E0U,0x9CE2AB57U,0x91A18D8EU,0x95609039U,
    0x8B27C03CU,0x8FE6DD8BU,0x82A5FB52U,0x8664E6E5U,
    0xBE2B5B58U,0xBAEA46EFU,0xB7A96036U,0xB3687D81U,
    0xAD2F2D84U,0xA9EE3033U,0xA4AD16EAU,0xA06C0B5DU,
    0xD4326D90U,0xD0F37027U,0xDDB056FEU,0xD9714B49U,
    0xC7361B4CU,0xC3F706FBU,0xCEB42022U,0xCA753D95U,
    0xF23A8028U,0xF6FB9D9FU,0xFBB8BB46U,0xFF79A6F1U,
    0xE13EF6F4U,0xE5FFEB43U,0xE8BCCD9AU,0xEC7DD02DU,
    0x34867077U,0x30476DC0U,0x3D044B19U,0x39C556AEU,
    0x278206ABU,0x23431B1CU,0x2E003DC5U,0x2AC12072U,
    0x128E9DCFU,0x164F8078U,0x1B0CA6A1U,0x1FCDBB16U,
    0x018AEB13U,0x054BF6A4U,0x0808D07DU,0x0CC9CDCAU,
    0x7897AB07U,0x7C56B6B0U,0x71159069U,0x75D48DDEU,
    0x6B93DDDBU,0x6F52C06CU,0x6211E6B5U,0x66D0FB02U,
    0x5E9F46BFU,0x5A5E5B08U,0x571D7DD1U,0x53DC6066U,
    0x4D9B3063U,0x495A2DD4U,0x44190B0DU,0x40D816BAU,
    0xACA5C697U,0xA864DB20U,0xA527FDF9U,0xA1E6E04EU,
    0xBFA1B04BU,0xBB60ADFCU,0xB6238B25U,0xB2E29692U,
    0x8AAD2B2FU,0x8E6C3698U,0x832F1041U,0x87EE0DF6U,
    0x99A95DF3U,0x9D684044U,0x902B669DU,0x94EA7B2AU,
    0xE0B41DE7U,0xE4750050U,0xE9362689U,0xEDF73B3EU,
    0xF3B06B3BU,0xF771768CU,0xFA325055U,0xFEF34DE2U,
    0xC6BCF05FU,0xC27DEDE8U,0xCF3ECB31U,0xCBFFD686U,
    0xD5B88683U,0xD1799B34U,0xDC3ABDEDU,0xD8FBA05AU,
    0x690CE0EEU,0x6DCDFD59U,0x608EDB80U,0x644FC637U,
    0x7A089632U,0x7EC98B85U,0x738AAD5CU,0x774BB0EBU,
    0x4F040D56U,0x4BC510E1U,0x46863638U,0x42472B8FU,
    0x5C007B8AU,0x58C1663DU,0x558240E4U,0x51435D53U,
    0x251D3B9EU,0x21DC2629U,0x2C9F00F0U,0x285E1D47U,
    0x36194D42U,0x32D850F5U,0x3F9B762CU,0x3B5A6B9BU,
    0x0315D626U,0x07D4CB91U,0x0A97ED48U,0x0E56F0FFU,
    0x1011A0FAU,0x14D0BD4DU,0x19939B94U,0x1D528623U,
    0xF12F560EU,0xF5EE4BB9U,0xF8AD6D60U,0xFC6C70D7U,
    0xE22B20D2U,0xE6EA3D65U,0xEBA91BBCU,0xEF68060BU,
    0xD727BBB6U,0xD3E6A601U,0xDEA580D8U,0xDA649D6FU,
    0xC423CD6AU,0xC0E2D0DDU,0xCDA1F604U,0xC960EBB3U,
    0xBD3E8D7EU,0xB9FF90C9U,0xB4BCB610U,0xB07DABA7U,
    0xAE3AFBA2U,0xAAFBE615U,0xA7B8C0CCU,0xA379DD7BU,
    0x9B3660C6U,0x9FF77D71U,0x92B45BA8U,0x9675461FU,
    0x8832161AU,0x8CF30BADU,0x81B02D74U,0x857130C3U,
    0x5D8A9099U,0x594B8D2EU,0x5408ABF7U,0x50C9B640U,
    0x4E8EE645U,0x4A4FFBF2U,0x470CDD2BU,0x43CDC09CU,
    0x7B827D21U,0x7F436096U,0x7200464FU,0x76C15BF8U,
    0x68860BFDU,0x6C47164AU,0x61043093U,0x65C52D24U,
    0x119B4BE9U,0x155A565EU,0x18197087U,0x1CD86D30U,
    0x029F3D35U,0x065E2082U,0x0B1D065BU,0x0FDC1BECU,
    0x3793A651U,0x3352BBE6U,0x3E119D3FU,0x3AD08088U,
    0x2497D08DU,0x2056CD3AU,0x2D15EBE3U,0x29D4F654U,
    0xC5A92679U,0xC1683BCEU,0xCC2B1D17U,0xC8EA00A0U,
    0xD6AD50A5U,0xD26C4D12U,0xDF2F6BCBU,0xDBEE767CU,
    0xE3A1CBC1U,0xE760D676U,0xEA23F0AFU,0xEEE2ED18U,
    0xF0A5BD1DU,0xF464A0AAU,0xF9278673U,0xFDE69BC4U,
    0x89B8FD09U,0x8D79E0BEU,0x803AC667U,0x84FBDBD0U,
    0x9ABC8BD5U,0x9E7D9662U,0x933EB0BBU,0x97FFAD0CU,
    0xAFB010B1U,0xAB710D06U,0xA6322BDFU,0xA2F33668U,
    0xBCB4666DU,0xB8757BDAU,0xB5365D03U,0xB1F740B4U};
   
  uint8_t data;
  uint32_t temporal,crc,total;
  int indice,ibyte,barrera;
  union{
    uint4 palabra;
    uint32_t vec[4];
  }sumas;

  total=MSC_V*NUMBETAS;
    
  for(barrera=0;barrera<3;barrera++)
    crc_thread_sums[barrera]=0xffffffffU;

  for(indice=0;indice<total;indice++){
    sumas.palabra=suma_global[indice];
    for(barrera=0;barrera<3;barrera++){
      temporal=sumas.vec[barrera];
      crc=crc_thread_sums[barrera];
      for(ibyte=0;ibyte<4;ibyte++){
	data=(uint8_t) (temporal&255U);
	temporal>>=8;
	crc = table_crc[data ^ ((crc >> 24) & 0xff)] ^ (crc << 8);
      }
      crc_thread_sums[barrera]=crc; 
    }
  }
      
}

void comprueba_resultados(long long int num_pasos,uint64_t semilla, uint32_t semilla_key)
{
  int indice,ib,site,barrera;
  double beta;
  
  double prob[3],error,empirical_prob,total;
  double esperado[3],dispersion[3];
  unsigned long long suma[3];
  uint32_t crc_thread_sums[3];
  union{
    uint4 palabra;
    uint32_t vec[4];
  }cuenta;
  FILE * Fout;

  if(NULL==(Fout=fopen("comprobacion_resultados_GPU.txt","wt")))
    print_and_exit("Problemas abriendo comprobacion_resultados_GPU.txt");
  
  fprintf(Fout,"Comprobando el generador en la GPU con semilla %llu\n", (unsigned long long) semilla);
  fprintf(Fout,"semilla_key=%u (binario:%u%u%u%u...)\n",semilla_key,(semilla_key>>31)&1U,
	 (semilla_key>>30)&1U,(semilla_key>>29)&1U,(semilla_key>>28)&1U);

  

  obten_crc(crc_thread_sums);

  indice=0;
  for(ib=0;ib<NUMBETAS;ib++){
    beta=betas[ib];
    for(barrera=0;barrera<3;barrera++){
      prob[barrera]=exp(-4.*beta*(barrera+1));
      esperado[barrera]=prob[barrera]*64.*((double) num_pasos);
      dispersion[barrera]=sqrt(prob[barrera]*(1.-prob[barrera])*64.*((double) num_pasos));
    }
    suma[0]=0; suma[1]=0;suma[2]=0;
    for(site=0;site<MSC_V;site++){
      cuenta.palabra=suma_global[indice];
      suma[0]+=suma_global[indice].x;
      suma[1]+=suma_global[indice].y;
      suma[2]+=suma_global[indice].z;
      for(barrera=0;barrera<3;barrera++){
	if(fabs((cuenta.vec[barrera]-esperado[barrera])/dispersion[barrera])>4.){
	  fprintf(Fout,"Gran desviacion ib=%d  site=%d barrera=%d: obtengo %u y espero %.14g (%.14g sigma)\n",
		 ib,site,4*(barrera+1),cuenta.vec[barrera],esperado[barrera],
		 (cuenta.vec[barrera]-esperado[barrera])/dispersion[barrera]);
	}
      }
      indice++;
    }
    total=((double) num_pasos*64.*MSC_V);

    for(barrera=0;barrera<3;barrera++){ 
      error=sqrt(prob[barrera]*(1-prob[barrera])/total);
      empirical_prob=((double) suma[barrera])/total;
      fprintf(Fout,"Densidad de bits en b_%d=%.14g (exact: %.14g, error=%.14g sigma) ib=%d\n",
	     4*(barrera+1),
	     empirical_prob,prob[barrera],(empirical_prob-prob[barrera])/error,
	     ib);
    }
  }
  fprintf(Fout,"Checksums:\n");
  for(barrera=0;barrera<3;barrera++){
    fprintf(Fout,"Barrera %d: suma_total=%llu,\t CRC-32 de las suma_global[tid]=%u\n",4*(barrera+1),
	   suma[barrera],crc_thread_sums[barrera]);
  }

  fclose(Fout);
}

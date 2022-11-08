/* File: func_dfuse.c */

/* ----------------------------------------------------------------- */
/* dfuse, a program to simulate hybridization and admixture.          */
/* Copyright (C) 2015, Doro Lindtke and Alex Buerkle                 */

/* This file is part of dfuse.                                       */

/* dfuse is free software: you can redistribute it and/or modify     */
/* it under the terms of the GNU General Public License as published */
/* by the Free Software Foundation, either version 3 of the License, */
/* or any later version.                                             */

/* dfuse is distributed in the hope that it will be useful,          */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of    */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the      */
/* GNU General Public License for more details.                      */

/* You should have received a copy of the GNU General Public License */
/* along with dfuse. If not, see <http://www.gnu.org/licenses/>.     */
/* ----------------------------------------------------------------- */

#include <string.h>
#include <ctype.h>
#include <limits.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_block.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_statistics_double.h>
#include "head_dfuse.h"



void usage(char * name){
  fprintf(stderr, "\nUsage: %s [options]\n", name);
  fprintf(stderr, "  -d  deme settings input filename\n");
  fprintf(stderr, "  -e  epistatic selection input filename\n");
  fprintf(stderr, "  -s  single-locus selection input filename\n");
  fprintf(stderr, "  -v  environment settings input filename\n");
  fprintf(stderr, "  -o  output filename [default = 'testing']\n");
  fprintf(stderr, "  -r  number of replicates [default = 1]\n");
  fprintf(stderr, "  -g  number of generations [default = 10]\n");
  fprintf(stderr, "  -G  print every x generations [default = 10]\n");
  fprintf(stderr, "  -O  print stats for selected loci every generation: 0 (no) or 1 (yes)\n");
  fprintf(stderr, "                      [default = 0]\n");
  fprintf(stderr, "  -m  migration rate [default = 0.05]\n");
  fprintf(stderr, "  -c  selection coefficient [default = 0.0]\n");
  fprintf(stderr, "  -a  selection against ancestral alleles on DMI ridge [default = 0.0]\n");
  fprintf(stderr, "  -R  recessivity for DMIs: 0 (dominant), 1 (co-dominant),\n");
  fprintf(stderr, "                      2 (partially recessive), 3 (fully recessive)\n");
  fprintf(stderr, "                      [default = 0]\n");
  fprintf(stderr, "  -S  selection stage 'viability', 'fertility', 'both' or\n");
  fprintf(stderr, "                      'none' [default = 'viability']\n");
  fprintf(stderr, "  -M  migration stage 'pollen', 'progeny' or 'both'\n");
  fprintf(stderr, "                      [default = 'progeny']\n");
  fprintf(stderr, "  -D  dispersal mode: 0 (infinite source) or 1 (finite source)\n");
  fprintf(stderr, "                      [default = 0]\n");
  fprintf(stderr, "  -l  marker loci per chromosome [default = 51]\n");
  /* fprintf(stderr, "  -f  fitness accumulation 'mult' or 'add' [default = 'mult']\n"); */
  fprintf(stderr, "  -Y  generation to start \"staining\" [default = 0 (no \"staining\")]\n");
  fprintf(stderr, "  -h  display options\n\n");
  exit(1);
}


/* ---------- */
/* read files */
/* ---------- */

/* deme file */
/* --------- */
/* first line, #demes, followed by line for adult capacities,
   line for progeny capacities */
void getdemes(char myfile[], int * ndemes, struct auxcont * aux){

  FILE *fp;
  gsl_vector_int * adultcap;
  gsl_vector_int * procap;
  int k;

  if ((fp = fopen(myfile,"r")) == NULL){
    fprintf(stderr, "Error: cannot open file %s\n", myfile);
    exit(1);
  }
  else {
    fscanf(fp, "%d", ndemes);
    adultcap = gsl_vector_int_calloc(*ndemes);
    procap = gsl_vector_int_calloc(*ndemes);

    gsl_vector_int_fscanf(fp, adultcap);
    gsl_vector_int_fscanf(fp, procap);

    if(dispersalmode==0){
      aux->cap_ad = gsl_vector_int_calloc(*ndemes);
      aux->cap_pr = gsl_vector_int_calloc(*ndemes);

      gsl_vector_int_memcpy(aux->cap_ad, adultcap);
      gsl_vector_int_memcpy(aux->cap_pr, procap);
    }
    else if(dispersalmode==1){ /* add peripheral demes */
      *ndemes += 2;
      aux->cap_ad = gsl_vector_int_calloc(*ndemes);
      aux->cap_pr = gsl_vector_int_calloc(*ndemes);

      gsl_vector_int_set(aux->cap_ad, 0, gsl_vector_int_get(adultcap, 0));
      gsl_vector_int_set(aux->cap_pr, 0, gsl_vector_int_get(procap, 0));

      for(k=0;k<(*ndemes)-2;k++){
	gsl_vector_int_set(aux->cap_ad, k+1, gsl_vector_int_get(adultcap, k));
	gsl_vector_int_set(aux->cap_pr, k+1, gsl_vector_int_get(procap, k));
      }

      gsl_vector_int_set(aux->cap_ad, (*ndemes)-1,
			 gsl_vector_int_get(adultcap, (*ndemes)-3));
      gsl_vector_int_set(aux->cap_pr, (*ndemes)-1,
			 gsl_vector_int_get(procap, (*ndemes)-3));
    }


    fclose(fp);
    gsl_vector_int_free(adultcap);
    gsl_vector_int_free(procap);
  }


  if(gsl_vector_int_max(aux->cap_pr) * kstrideInd > INT_MAX){
    fprintf(stderr,"Error: progeny capacity of %d can't be handled\n",
	    gsl_vector_int_max(aux->cap_pr));
    exit(1);
  }
  if(gsl_vector_int_max(aux->cap_ad) * kstrideInd > INT_MAX){
    fprintf(stderr,"Error: adult capacity of %d can't be handled\n",
	    gsl_vector_int_max(aux->cap_ad));
    exit(1);
  }
  if(gsl_vector_int_max(aux->cap_ad) * MeanFecund > INT_MAX){
    fprintf(stderr,"Error: adult capacity of %d and fecundity of %.1f can't be handled\n",
	    gsl_vector_int_max(aux->cap_ad), (double) MeanFecund * 2);
    exit(1);
  }
  fprintf(stderr, "Reading demefile %s done\n", myfile);
}


/* epistatic selection file */
/* ------------------------ */
void getepi(char myfile[], int * nepiloci,
	    struct epiloci * epi, const double sel, const double ridge,
	    const int recessivity){


  int ncompl, nplets;
  gsl_vector_int * type;
  gsl_vector_int * ancallele;

  int i;
  int currcomp = 0, currcompcount = 0;
  FILE *fp;
  if ((fp = fopen(myfile,"r")) == NULL){
    fprintf(stderr, "Error: cannot open file %s\n", myfile);
    exit(1);
  }
  else {
    fscanf(fp, "%d%d%d", nepiloci, &ncompl, &nplets);

    if (*nepiloci != (ncompl * nplets)) {
      fprintf(stderr,
	      "Error: %d selected loci do not match %d complexes with %d-plets\n",
	      *nepiloci, ncompl, nplets);
      exit(1);
    }
    else if (nplets < 2) {
      fprintf(stderr,
	      "Error: cannot calculate epistatic interactions with %d-plets\n",
	      nplets);
      exit(1);
    }
    else {
      /* fprintf(stderr, "%d epistatic loci, %d complexes, %d-plets\n",  */
      /* 	      *nepiloci, ncompl, nplets); */
      ;
    }

    epi->ncomplexes = ncompl;
    epi->xplets = nplets;

    epi->locuscount = gsl_vector_int_calloc(*nepiloci);
    epi->complex = gsl_vector_int_calloc(*nepiloci);
    epi->compcount = gsl_vector_int_calloc(*nepiloci);
    epi->type = gsl_vector_char_calloc(*nepiloci);
    epi->chr = gsl_vector_int_calloc(*nepiloci);
    epi->posdecimal = gsl_vector_calloc(*nepiloci);
    epi->ancallele = gsl_vector_char_calloc(*nepiloci);

    type = gsl_vector_int_calloc(*nepiloci);
    ancallele = gsl_vector_int_calloc(*nepiloci);


    gsl_vector_int_fscanf(fp, epi->locuscount);
    gsl_vector_int_fscanf(fp, epi->complex);
    gsl_vector_int_fscanf(fp, epi->compcount);
    gsl_vector_int_fscanf(fp, type);
    gsl_vector_int_fscanf(fp, epi->chr);
    gsl_vector_fscanf(fp, epi->posdecimal);
    gsl_vector_int_fscanf(fp, ancallele);


    for (i=0; i < *nepiloci; i++){
      if (gsl_vector_int_get(type,i)==0){
	gsl_vector_char_set(epi->type,i,'d');
      }
      else if (gsl_vector_int_get(type,i)==1){
	gsl_vector_char_set(epi->type,i,'p');
      }
      else {
	gsl_vector_char_set(epi->type,i,'w');
	fprintf(stderr, "Warning: epistatic selection type unknown\n");
      }
    }

    for (i=0; i < *nepiloci; i++){
      if (gsl_vector_int_get(ancallele,i)==0){
	gsl_vector_char_set(epi->ancallele,i,'a');
      }
      else if (gsl_vector_int_get(ancallele,i)==1){
	gsl_vector_char_set(epi->ancallele,i,'b');
      }
      else {
	gsl_vector_char_set(epi->ancallele,i,'w');
	fprintf(stderr, "Warning: ancestral allele for epistasis unknown\n");
      }
    }


    fclose(fp);

    /* some checks... */
    for (i=0; i < *nepiloci; i++){
      /* check that locuscount, complexnumber and complexcount in input file are correct */
      if (gsl_vector_int_get(epi->locuscount, i) != i){
	fprintf(stderr, "Error: locus counter in %s is wrong\n",myfile);
	exit(1);
      }
      if (gsl_vector_int_get(epi->complex, i) != currcomp){
	fprintf(stderr, "Error: complex number in %s is wrong\n",myfile);
	exit(1);
      }
      if (gsl_vector_int_get(epi->compcount, i) == currcompcount){
	currcompcount++;
      }
      else {
	fprintf(stderr, "Error: complex counter in %s is wrong\n",myfile);
	exit(1);
      }
      if ((i + 1) % nplets == 0){
	currcomp++;
	currcompcount = 0;
      }
      /* error if position > 1 or < 0 */
      if (gsl_vector_get(epi->posdecimal, i) < 0 || gsl_vector_get(epi->posdecimal, i) > 1){
	fprintf(stderr, "Error: locus position in %s not between 0.0 and 1.0 Morgans\n",myfile);
	exit(1);
      }
      /* error if chromosome > (khomologsN - 1) or < 0 (when starting at 0)*/
      if (gsl_vector_int_get(epi->chr, i) < 0 || gsl_vector_int_get(epi->chr, i) >
	  (khomologsN - 1)){
	fprintf(stderr, "Error: chromosome number in %s impossible\n",myfile);
	exit(1);
      }
    }

  }
  /* set fitness matrices for epistasis */
  setfitmat(sel, epi->dmimat, epi->pathmat, ridge, recessivity);

  gsl_vector_int_free(type);
  gsl_vector_int_free(ancallele);

  fprintf(stderr, "Reading epifile %s done\n", myfile);

}

/* single-locus selection file */
/* --------------------------- */
int getsel(char myfile[], int * nselloci, struct selloci * sel){

  gsl_vector_int * type;
  int i, env = 0;
  FILE *fp;

  if ((fp = fopen(myfile,"r")) == NULL){
    fprintf(stderr, "Error: cannot open file %s\n", myfile);
    exit(1);
  }
  else {
    fscanf(fp, "%d", nselloci);

    /* fprintf(stderr, "%d selected loci\n", *nselloci); */

    sel->locuscount = gsl_vector_int_calloc(*nselloci);
    sel->type = gsl_vector_char_calloc(*nselloci);
    sel->chr = gsl_vector_int_calloc(*nselloci);
    sel->posdecimal = gsl_vector_calloc(*nselloci);
    sel->valaa = gsl_vector_calloc(*nselloci);
    sel->valab = gsl_vector_calloc(*nselloci);
    sel->valbb = gsl_vector_calloc(*nselloci);

    type = gsl_vector_int_calloc(*nselloci);


    gsl_vector_int_fscanf(fp, sel->locuscount);
    gsl_vector_int_fscanf(fp, type);
    gsl_vector_int_fscanf(fp, sel->chr);
    gsl_vector_fscanf(fp, sel->posdecimal);
    gsl_vector_fscanf(fp, sel->valaa);
    gsl_vector_fscanf(fp, sel->valab);
    gsl_vector_fscanf(fp, sel->valbb);

    for (i=0; i < *nselloci; i++){
      if (gsl_vector_int_get(type,i)==1){
	gsl_vector_char_set(sel->type,i,'e');
	env++; /* count environmental selected loci */
      }
      else if (gsl_vector_int_get(type,i)==2){
	gsl_vector_char_set(sel->type,i,'o');
      }
      else if (gsl_vector_int_get(type,i)==3){
	gsl_vector_char_set(sel->type,i,'u');
      }
      else {
	gsl_vector_char_set(sel->type,i,'w');
	fprintf(stderr, "Warning: selection type in %s unknown\n",myfile);
      }
    }

   fclose(fp);

   /* some checks... */
   for (i=0; i < *nselloci; i++){
     /* check that locuscount in input file is correct */
     if (gsl_vector_int_get(sel->locuscount, i) != i){
       fprintf(stderr, "Error: locus counter in %s is wrong\n",myfile);
       exit(1);
     }
     /* error if position > 1 or < 0 */
     if (gsl_vector_get(sel->posdecimal, i) < 0 || gsl_vector_get(sel->posdecimal, i) > 1){
       fprintf(stderr, "Error: locus position in %s not between 0.0 and 1.0 Morgans\n",myfile);
       exit(1);
     }
     /* error if chromosome > (khomologsN - 1) or < 0 */
     if (gsl_vector_int_get(sel->chr, i) < 0 || gsl_vector_int_get(sel->chr, i) >
	 (khomologsN - 1)){
       fprintf(stderr, "Error: chromosome number in %s impossible\n",myfile);
       exit(1);
     }
     /* error if genotype values  > 1 or < 0 */
     if (gsl_vector_get(sel->valaa, i) < 0 || gsl_vector_get(sel->valaa, i) > 1){
       fprintf(stderr, "Error: value for genotype 'aa' not between 0.0 and 1.0\n");
       exit(1);
     }
     if (gsl_vector_get(sel->valab, i) < 0 || gsl_vector_get(sel->valab, i) > 1){
       fprintf(stderr, "Error: value for genotype 'ab' not between 0.0 and 1.0\n");
       exit(1);
     }
     if (gsl_vector_get(sel->valbb, i) < 0 || gsl_vector_get(sel->valbb, i) > 1){
       fprintf(stderr, "Error: value for genotype 'bb' not between 0.0 and 1.0\n");
       exit(1);
     }
   }

  }
  gsl_vector_int_free(type);

  fprintf(stderr, "Reading selfile %s done\n", myfile);

  return(env);
}



/* read enviloci x deme file */
/* ------------------------- */
void getenvi(char myfile[], const int ndemes, const int nenviloci, struct auxcont * aux){

  int i, k, loci, patches;
  FILE *fp;
  gsl_matrix * envimat;

  if ((fp = fopen(myfile,"r")) == NULL){
    fprintf(stderr, "Error: cannot open file %s\n", myfile);
    exit(1);
  }
  else {
    fscanf(fp, "%d%d", &loci, &patches);
    if(loci != nenviloci){
      fprintf(stderr, "Error: conflicting number of loci in %s\n",myfile);
      exit(1);
    }
    if(dispersalmode==0 && patches != ndemes){
      fprintf(stderr, "Error: conflicting number of demes in %s\n",myfile);
      exit(1);
    }
    if(dispersalmode==1 && patches+2 != ndemes){
      fprintf(stderr, "Error: conflicting number of demes in %s\n",myfile);
      exit(1);
    }

    aux->envi = gsl_matrix_calloc(loci * 2, ndemes);
    envimat = gsl_matrix_calloc(loci * 2, patches);

    if (gsl_matrix_fscanf(fp, envimat) != 0 ){
      fprintf(stderr, "Error with scanning environmental parameters\n");
      exit(1);
    }

    fclose(fp);

    /* some checks... */
    for (i=0;i<(loci*2);i++){
      for (k=0;k<patches;k++){
	if (gsl_matrix_get(envimat, i, k) < 0 || gsl_matrix_get(envimat, i, k) > 1){
	  fprintf(stderr, "Error: environmental parameters not between 0.0 and 1.0\n");
	  exit(1);
	}
      }
    }

  }

  /* copy temp matrix to real one */
  if(dispersalmode==0){
    gsl_matrix_memcpy(aux->envi,envimat);
  }
  else if (dispersalmode==1){
    for(i=0;i<(loci*2);i++){
      gsl_matrix_set(aux->envi, i,0, gsl_matrix_get(envimat, i,0));
      for(k=0;k<patches;k++){
	gsl_matrix_set(aux->envi,i, k+1, gsl_matrix_get(envimat,i, k));
      }
      gsl_matrix_set(aux->envi, i,patches+1, gsl_matrix_get(envimat, i,patches-1));
    }
  }

  gsl_matrix_free(envimat);

  fprintf(stderr, "Reading envifile %s done\n", myfile);
}

/* ----------------- */
/* setting things up */
/* ----------------- */

void defaultenvi(const double selstrength,
		 const int patches, const int loci, struct auxcont * aux){

  int i, k;
  gsl_vector * habitat;
  gsl_vector * sel;

  aux->envi = gsl_matrix_calloc(loci * 2, patches);
  habitat = gsl_vector_alloc(patches);
  sel = gsl_vector_alloc(patches);

  gsl_vector_set_all(sel,selstrength);
  if (patches == 1){ /* if only one deme, hybrid habitat will be assigned */
    gsl_vector_set(habitat, 0, 0.5);
  }
  else {
    if(dispersalmode==0){
      for (k=0;k<patches;k++){ /* gradual change from 0 to 1 */
	gsl_vector_set(habitat, k, (double) k / (double) (patches - 1));
      }
    }
    else if(dispersalmode==1){ /* peripheral demes 0 or 1, middle demes gradual change */
      gsl_vector_set(habitat, 0, 0.0);
      for (k=0;k<patches-2;k++){ /* gradual change from 0 to 1 */
	gsl_vector_set(habitat, k+1, (double) k / (double) (patches - 3));
      }
      gsl_vector_set(habitat, patches-1, 1.0);
    }
  }

  for (i=0;i<(loci*2);i++){
    if(i % 2 == 0){ /* set environment */
      gsl_matrix_set_row(aux->envi, i, habitat);
    }
    else if(i % 2 == 1){ /* set sel coeff */
      gsl_matrix_set_row(aux->envi, i, sel);
    }
  }

  gsl_vector_free(habitat);
  gsl_vector_free(sel);
}




/* set migration matrix (stepping stone as default) */
void setmigmat(const int ndemes, struct auxcont * aux, const double mig){

  int i, j;

  /* rows for recipient deme, cols for source deme (+2 for pure species
     at each side of the deme-chain)  */
  aux->migmat = gsl_matrix_calloc(ndemes, ndemes + 2);

  if (dispersalmode==0){ /* infinite source */
    for (i=0; i< ndemes; i++){ /* recipient */
      for (j=0; j< (ndemes + 2); j++){ /* source */
	if (i+1 == j){
	  gsl_matrix_set(aux->migmat,i,j,(1.0 - mig)); /* residents */
	}
	else if (i == j || i+2 == j){
	  gsl_matrix_set(aux->migmat,i,j,mig/2); /* migrants */
	}
      }
    }
  }

  else if (dispersalmode==1){ /* finite source; no migration from or to source demes */
    for (i=0; i< ndemes; i++){ /* recipient */
      for (j=0; j< (ndemes + 2); j++){ /* source */
	if (j==0 || j==(ndemes + 1)){
	  ;
	}
	else if (i+1 == j && (j==1 || j== (ndemes))){
	  gsl_matrix_set(aux->migmat,i,j,(1.0 - (mig/2))); /* residents at peripheres*/
	}
	else if (i+1 == j && (j>1 || j< (ndemes))){
	  gsl_matrix_set(aux->migmat,i,j,(1.0 - mig)); /* residents in middle */
	}
	else if (i == j || i+2 == j){
	  gsl_matrix_set(aux->migmat,i,j,mig/2); /* migrants */
	}
      }
    }
  }

}

/* set fitness matrix */
void setfitmat(const double sel, double dmi[][3], double path[][3], const double ridge,
	       const int recessivity){

  int i, j;
  double h0, h1, h2;
  h2 = sel;

  if (recessivity==1){ /* co-dominant DMIs */
    h0 = sel/4;
    h1 = sel/2;
  }
  else if (recessivity==2){ /* partially recessive DMIs */
    h0 = 0.0;
    h1 = sel/2;
  }
  else if (recessivity==3){ /* recessive DMIs */
    h0 = 0.0;
    h1 = 0.0;
  }
  else { /* dominant DMIs */
    h0 = sel;
    h1 = sel;
  }

  /* fitness matrix to specify fitness according to number of derived alleles
     in pairwise calculations */
  double mydmi[3][3] = {
    {1.0 - ridge,    1.0 - ridge/2,  1.0      },
    {1.0 - ridge/2,  1.0 - h0,      1.0 - h1},
    {1.0,            1.0 - h1,      1.0 - h2},
  };


  /* fitness matrix to specify fitness according to number of mismatching alleles
     in pairwise calculations for linear pathways */
  double mypath[3][3] = {
    {1.0,         1.0 - sel/2,   1.0 - sel  },
    {1.0 - sel/2, 1.0,           1.0 - sel/2},
    {1.0 - sel,   1.0 - sel/2,   1.0        },
  };

  for(i=0;i<3;i++){
    for(j=0;j<3;j++){
      dmi[i][j] = mydmi[i][j];
      path[i][j] = mypath[i][j];
    }
  }

  /* printf("\nDMI fitness matrix:\n"); */
  /* for(i=0; i<3; i++){ */
  /*   for (j=0; j< 3; j++){ */
  /*     printf("%7.3f  ", dmi[i][j]); */
  /*   } */
  /*   printf("\n"); */
  /* } */

}



/* set deme grid */
/* first struct is array of struct, second struct is pointer to one struct */
void init_deme(struct deme * grid, struct auxcont * aux, const int nenviloci){

  int i, j;
  unsigned int popbreak;

  /* allocate memory and set some vals */
  for (i=0; i<ndemes; i++){
    grid[i].cap_ad = gsl_vector_int_get(aux->cap_ad, i);
    if (grid[i].cap_ad < 1){
      fprintf(stderr,"Error: adult capacity cannot be smaller than 1\n");
      exit(1);
    }
    grid[i].cap_pr = gsl_vector_int_get(aux->cap_pr, i);
    /* capacity * kstrideInd = all chromosomes in deme */
    grid[i].gChromoNpop = gsl_vector_int_get(aux->cap_ad, i) * kstrideInd;
    grid[i].N = grid[i].cap_ad; /* pop size */
    grid[i].Npr = 0; /* number of progeny */
    grid[i].fitness = gsl_vector_alloc(grid[i].cap_ad);
    gsl_vector_set_all(grid[i].fitness, 0.0);
    grid[i].stats = gsl_matrix_calloc(grid[i].cap_ad, 2); /* for het and hi */
    grid[i].junctions = gsl_vector_alloc(grid[i].cap_ad); /* for # junctions per ind */
    grid[i].sub_cur_generation = malloc((gsl_vector_int_get(aux->cap_ad, i) *
					 kstrideInd) * sizeof(struct chromosome));
    if(grid[i].sub_cur_generation == NULL){
      printf("Error in allocating memory for deme %d\n", i);
    }
    grid[i].dummyHead_sub_next_gen = malloc(sizeof(struct dyn_chromosome));

    /* environment, if specified */
    if(nenviloci > 0){
      grid[i].environment = gsl_vector_alloc(nenviloci);
      grid[i].envselstrength = gsl_vector_alloc(nenviloci);
      for (j=0;j<nenviloci;j++){
	gsl_vector_set(grid[i].environment, j, gsl_matrix_get(aux->envi, j*2, i));
	gsl_vector_set(grid[i].envselstrength, j, gsl_matrix_get(aux->envi, j*2 + 1, i));
      }
    }

  }

  /* initialize chromosomes in demes */
  /* for now, fill parentals from both sides to middle of metapop grid */
  for (i = 0; i < ndemes/2; i++){
    popbreak = grid[i].gChromoNpop;
    /* fill parental species (all type 'a') */
    pop_init(grid[i].sub_cur_generation, popbreak, grid[i].gChromoNpop);
  }
  if (ndemes % 2 == 0){
    for (i = ndemes-1 ; i >= ndemes/2; i--){
      popbreak = 0;
      /* fill parental species (all type 'b') */
      pop_init(grid[i].sub_cur_generation, popbreak, grid[i].gChromoNpop);
    }
  }
  else if (ndemes % 2 == 1){
    for (i = ndemes-1 ; i >= ndemes/2; i--){
      popbreak = (i==ndemes/2) ? (int)(grid[i].cap_ad/2) * kstrideInd : 0;
      pop_init(grid[i].sub_cur_generation, popbreak, grid[i].gChromoNpop);
    }
  }

}


/* sets array of individuals/homologs to starting conditions/genotypes */
void pop_init(struct chromosome * cur_generation, unsigned int popbreak,
	      unsigned int gChromoNpop){
  unsigned int k;

  /* Parental species 'a', left of grid */
  for (k=0; k < popbreak; k++){
    init_chromosome (&cur_generation[k], 1); /* allocate memory for chromosome k, 1 block */
    gsl_vector_char_set(cur_generation[k].allele, 0,'a'); /* set element 0 to 'a' */
    gsl_vector_float_set(cur_generation[k].mapposition, 0, 1); /* set element 0 to 1.000 */
    cur_generation[k].nblocks = 1;
  }

  /* Parental species 'b', right of grid */
  for (k=popbreak; k < gChromoNpop; k++){
    init_chromosome (&cur_generation[k], 1);
    gsl_vector_char_set(cur_generation[k].allele, 0,'b');
    gsl_vector_float_set(cur_generation[k].mapposition, 0, 1);
    cur_generation[k].nblocks = 1;
  }
}


/* --------------------- */
/* simulation essentials */
/* --------------------- */

int init_chromosome (struct chromosome * chromo, int nblocks){
  chromo->allele = gsl_vector_char_alloc(nblocks);
  chromo->mapposition = gsl_vector_float_alloc(nblocks);
  return(0);
}


int init_dyn_chromosome (struct dyn_chromosome * chromo, int nblocks){
  chromo->allele = gsl_vector_char_alloc(nblocks);
  chromo->mapposition = gsl_vector_float_alloc(nblocks);
  return(0);
}

int find_position (double position, const struct chromosome * chromo){
/* finds the index for the block containing the locus position */
  int i;
  for(i=0; i< chromo->nblocks; i++){
    if(position <= gsl_vector_float_get(chromo->mapposition, i) ){
      break;
    }
  }
  return(i);
}



void recombine (const struct chromosome * first,
		const struct chromosome * second,
		struct chromosome * gamete){

  int i, k, tmpnblocks, run, elementCount, end, ncrossovers;

  gsl_vector * chiasmaLoc;
  gsl_vector_char * tmpAllele;
  gsl_vector_char_view subsetAllele;
  gsl_vector_float * tmpMapposition;
  gsl_vector_float_view subsetMapposition;
  gsl_vector_int * firstindex, * secondindex;
  /* ASSUMPTIONS ---- always start gamete on first homolog */

  /* parameter of Poisson is a compile-time variable */
  ncrossovers = gsl_ran_poisson(r, kchiasmaN);

  /* either no crossovers occur, or some (given by possibly modified
     ncrossovers) do occur */
  if(ncrossovers == 0){
    /* in this case one or the other homolog should be given to progeny */
    if(gsl_rng_uniform(r) < 0.5){
      init_chromosome(gamete, first->nblocks); /* allocate memory */
      copy_homolog(gamete, first);
    }
    else{
      init_chromosome(gamete, second->nblocks);
      copy_homolog(gamete, second);
    }
  }
  else{
    chiasmaLoc = gsl_vector_alloc(ncrossovers);
    /* figure out how many crossovers actually are going to occur, and
       place them */

    for(i=0;i<ncrossovers;i++){
      gsl_vector_set(chiasmaLoc, i, gsl_rng_uniform(r) );
    }
    gsl_sort_vector(chiasmaLoc);

    tmpnblocks = GSL_MAX_INT(first->nblocks, second->nblocks) * (ncrossovers + 1);
    /* note this can be 1*(1+1), so the plus one is necessary
       (or when there is little overlap among blocks from the two homologs) */
    tmpAllele = gsl_vector_char_alloc(tmpnblocks);
    tmpMapposition = gsl_vector_float_alloc(tmpnblocks);

    firstindex = gsl_vector_int_alloc(ncrossovers);
    secondindex = gsl_vector_int_alloc(ncrossovers);

    for(i=0;i<ncrossovers;i++){
      gsl_vector_int_set(firstindex, i,
			 find_position( gsl_vector_get(chiasmaLoc, i), first) );
      /* find_position returns ID of last block before new breakpoint
	 (-> block where chiasma occurs) */
      gsl_vector_int_set(secondindex, i,
			 find_position( gsl_vector_get(chiasmaLoc, i), second));
    }

    /* build gamete */
    run = elementCount = end = 0;

    /*  fill run 0, with check to see whether chiasma occurs in first block */
    if(0 != gsl_vector_int_get(firstindex, 0)){
      /* first chiasma is not in first block */
      end = gsl_vector_int_get(firstindex, 0); /* block until chromosome can simply be copied */

      fillgamete(tmpAllele, tmpMapposition, first, 0,
		 end, &elementCount); /* copy first part of chromosome into tmp */
      /* append block with chiasma (copy next block until point of chisma) */
      appendChiasma(tmpAllele, elementCount, gsl_vector_char_get(first->allele, end),
		    tmpMapposition, gsl_vector_get(chiasmaLoc, run));
      elementCount++;

    }
    else{ /* if first chiasma is in first block */
      appendChiasma(tmpAllele, elementCount, gsl_vector_char_get(first->allele, 0),
		     tmpMapposition, gsl_vector_get(chiasmaLoc, run));
      elementCount++;
    }
    run++; /* go to next chiasma */

    for(; run<ncrossovers; run++){
      if(GSL_IS_ODD(run)){ /* run is odd -> copy from second homolog */
	/* block until chromosome can simply be copied */
	end = gsl_vector_int_get(secondindex, run);

	fillgamete(tmpAllele, tmpMapposition, second,
		   gsl_vector_int_get(secondindex, run-1), end, &elementCount);
	/* append block with chiasma */
	appendChiasma(tmpAllele, elementCount,
		       gsl_vector_char_get(second->allele, end),
		       tmpMapposition, gsl_vector_get(chiasmaLoc, run));
	elementCount++;

      }
      else { /* run is even -> copy from first homolog */
	end = gsl_vector_int_get(firstindex, run);
	fillgamete(tmpAllele, tmpMapposition, first,
		   gsl_vector_int_get(firstindex, run-1), end, &elementCount);
	/* append block with chiasma */
	appendChiasma(tmpAllele, elementCount, gsl_vector_char_get(first->allele, end),
		      tmpMapposition, gsl_vector_get(chiasmaLoc, run));
	elementCount++;
      }
    }

    /* fill final run (copy remaining part of chromosome) */
    if(GSL_IS_ODD(run)){
      fillgamete(tmpAllele, tmpMapposition, second,
		 gsl_vector_int_get(secondindex, run-1),
		 second->nblocks, &elementCount);
    }
    else {
      fillgamete(tmpAllele, tmpMapposition, first,
		 gsl_vector_int_get(firstindex, run-1),
		 first->nblocks, &elementCount);
    }

    /* simplify both tmp vectors in parallel, verified that this works */
    /* this will remove junctions that are between blocks of identical ancestry */
    tmpnblocks = elementCount;

    for(k=0, i=1; i<elementCount; i++){
      if(gsl_vector_char_get(tmpAllele, i) == gsl_vector_char_get(tmpAllele, i-1)){
	k++; /* vector offset to move items into correct position*/
	tmpnblocks--;
      }
      gsl_vector_char_set(tmpAllele, i-k, gsl_vector_char_get(tmpAllele, i) );
      gsl_vector_float_set(tmpMapposition, i-k, gsl_vector_float_get(tmpMapposition, i) );
    }

    /*  and store the simplified data in the persistant gamete variable */
    init_chromosome(gamete, tmpnblocks); /* allocate memory */

    subsetAllele = gsl_vector_char_subvector(tmpAllele, 0, tmpnblocks);
    subsetMapposition = gsl_vector_float_subvector(tmpMapposition, 0, tmpnblocks);
    gsl_vector_char_memcpy(gamete->allele, &subsetAllele.vector );
    gsl_vector_float_memcpy(gamete->mapposition, &subsetMapposition.vector );

    gamete->nblocks = tmpnblocks;

    gsl_vector_char_free(tmpAllele);
    gsl_vector_float_free(tmpMapposition);
    gsl_vector_int_free(firstindex);
    gsl_vector_int_free(secondindex);
    gsl_vector_free(chiasmaLoc);
  }
}

int fillgamete (gsl_vector_char * allele, gsl_vector_float * position,
		const struct chromosome * homolog,
		const int begin, const int end, int * elementCount){
  /* copies parts of the chromosome */
  int i;

  for (i=begin; i<end; i++){
    gsl_vector_char_set(allele, *elementCount,
			gsl_vector_char_get(homolog->allele, i ));
    gsl_vector_float_set(position, *elementCount,
		   gsl_vector_float_get(homolog->mapposition, i ));
    *elementCount = * elementCount + 1;
  }
  return(0);
}

int appendChiasma (gsl_vector_char * allele,
		   const int elementCount,
		   char alleletype,
		   gsl_vector_float * mapposition,
		   double location ){
  gsl_vector_char_set(allele, elementCount, alleletype  );
  gsl_vector_float_set(mapposition, elementCount, location );

  return(0);
}

int storehomolog (struct dyn_chromosome * homolog, struct chromosome * gamete){
  gsl_vector_char_memcpy(homolog->allele, gamete->allele);
  gsl_vector_float_memcpy(homolog->mapposition, gamete->mapposition);
  homolog->nblocks = gamete->nblocks;
  return(0);
}

int storehomologViab (struct chromosome * homolog, struct chromosome * gamete){
  gsl_vector_char_memcpy(homolog->allele, gamete->allele);
  gsl_vector_float_memcpy(homolog->mapposition, gamete->mapposition);
  homolog->nblocks = gamete->nblocks;
  return(0);
}


/* ---------- */
/* life cycle */
/* ---------- */

/* females (each resident in deme can be female) */
/* assign pot. #progeny from poisson for each female
   (fertility will affect poisson parameter) */
/* for each mom, store ID times #progeny in vector for potential progeny;
   shuffle vector */
/* make maternal gamete for each potential progeny, and... */

/* ...choose dad, migrant^ or resident (no selfing) (-> pollen dispersal) */
/* test dad's fecundity, make gamete */

/* make zygote, test zygote (relative fitness, viability) */
/* add zygote to dyn_chromosome */
/* while adding, also migrate (-> seed or progeny dispersal) */
/* end if no more pot. #progeny, or progeny carrying capacity filled */

/* copy to cur_generation later with copynext2cur(),
   reduce to adult capacity (random or competition ->
   this would be double-selected) */

/* ^ for migrants: marginal demes will receive migrants from pure species */

/* --------------------------------------------------------------------- */

void mating (struct deme * grid,
	     struct auxcont * aux,
	     char * selstage, char * migstage,
	     const int nallselloci,
	     const struct epiloci * epi, const struct selloci * sel,
	     const int nepiloci, const int nselloci, const double selstrength){

  int k, m, source, sink;
  unsigned int tmpnpro;
  unsigned int mat_index, pat_index, tmp, nprocap, nproposs, homolog, i, j;
  int found; /* 1 indicates parent is found */
  struct chromosome paternalgamete, maternalgamete;
  struct chromosome potentialzygote[(2*kstrideInd)];
  gsl_vector * migration;
  gsl_vector * migration_copy;
  gsl_vector_uint * sourcevec;
  gsl_vector_uint * sinkvec;
  gsl_vector_uint * procntr; /* progeny counter for each female */
  gsl_vector_uint * gameteID; /* ids of females for potential gametes */
  gsl_vector_int * geno; /* to store genotypes (0,1,2 'a' alleles) at selected sites */
  double zygotefitness, p;
  char allele = 'c';

  migration = gsl_vector_calloc(ndemes + 2);
  migration_copy = gsl_vector_calloc(ndemes + 2);
  sourcevec = gsl_vector_uint_calloc(ndemes + 2);
  sinkvec = gsl_vector_uint_calloc(ndemes + 2);
  if(nallselloci > 0){
    geno = gsl_vector_int_alloc(nallselloci);
  }

  /* for deme */
  for (k=0; k<ndemes; k++){
    /* migration; specify source deme of parent */
    gsl_matrix_get_row(migration, aux->migmat, k);

    /* set up reproductive females */
    /* all inds in deme act as females (no migrants) */
    /* --------------------------------------------- */
    if (grid[k].N > 0) {
      procntr = gsl_vector_uint_calloc(grid[k].N);
      nprocap = grid[k].cap_pr; /* max possible progeny in deme given capacity */
      nproposs = 0; /* max possible progeny in deme given sum of female fecundity */

      /* take females, assign number of progeny from poisson distr,
	 poisson param multiplied by female fecundity (fitness) */
      for (i=0; i<grid[k].N; i++){
	if (!strcmp(selstage,"viability") || !strcmp(selstage,"none")){ /* no fertility selection */
	  tmpnpro = gsl_ran_poisson(r, (double) MeanFecund);
	}
	else { /* with fertility selection */
	  tmpnpro = gsl_ran_poisson(r, (double) MeanFecund * gsl_vector_get(grid[k].fitness, i));
	}
	gsl_vector_uint_set(procntr, i, tmpnpro);
	nproposs += tmpnpro;
      }
      /* store all mat_index for #progeny in vector and move through them later */
      gameteID = gsl_vector_uint_calloc(nproposs);
      tmpnpro = 0;
      for (i=0; i<grid[k].N; i++){
	for (j=0; j < gsl_vector_uint_get(procntr, i); j++){
	  gsl_vector_uint_set(gameteID, tmpnpro, i);
	  tmpnpro++;
	}
      }
      gsl_ran_shuffle(r, gameteID->data, nproposs, sizeof(unsigned int));

      gsl_vector_uint_free(procntr);

      /* pop size 1, no pollen migration, no selfing -> don't make any progeny */
      /* (immigration from pure species will still be possible) */
      if((!strcmp(migstage, "progeny") || gsl_vector_get(migration, k+1)==1.0) && grid[k].N == 1){
	nproposs = 0;
      }

      /* make zygotes */
      /* ------------ */
      while(nprocap>0 && nproposs>0){

	/* take mom, make maternal gamete below */
	nproposs--; /* always decrement if maternal gamete made */
	mat_index = gsl_vector_uint_get(gameteID, nproposs); /* starts at end, but does not matter */

	/* choose dad, make gamete below */
	/* dad can be migrant */
	/* ------------------ */
	found = 0;
	source = -1;

	if (!strcmp(migstage, "progeny") || gsl_vector_get(migration, k+1)==1.0){
	  /* no pollen migration */
	  source = k+1;
	}
	else { /* pollen migration */
	  /* this assumes that pollen is always unlimited (no matter what deme's N is) */
	  gsl_vector_memcpy(migration_copy, migration);
	  if(grid[k].N == 1){ /* no selfing: don't sample pollen from same deme */
	    gsl_vector_set(migration_copy, k+1, 0.0); /* remaining probs will automatically
							 be normalized below */
	  }

	  gsl_ran_multinomial(r, ndemes+2, 1, migration_copy->data, sourcevec->data);

	  /* source will have zeros and one 1 */
	  for (m=0; m<ndemes+2; m++){
	    if (gsl_vector_uint_get(sourcevec,m) == 1){
	      source = m;
	    }
	  }
	}


	if (source == -1){
	  fprintf(stderr, "Error in assigning source deme of paternal gamete\n");
	  exit(1);
	}
	else if (source == 0){
	  /* gamete from parental 0 */
	  found = 1; /* full fertility */
	}
	else if (source == ndemes+1){
	  /* gamete from parental 1 */
	  found = 1; /* full fertility */
	}
	else if (source - 1 == k){
	  /* check that no selfing (only possible if male is resident) */
	  while(found == 0){
	    tmp = gsl_rng_uniform_int(r, grid[k].N);

	    if(tmp != mat_index){ /* disallow selfing */
	      if(!strcmp(selstage, "viability") || !strcmp(selstage, "none")){
		pat_index = tmp; /* no fertility selection */
		found = 1;
	      }
	      else if(gsl_rng_uniform(r) < gsl_vector_get(grid[k].fitness, tmp) ){
		pat_index = tmp;
		found = 1;
	      }
	    }
	  }
	}
	else {
	  /* migrant; fertility selection applies in source deme */
	  while(found == 0){
	    tmp = gsl_rng_uniform_int(r, grid[source - 1].N);

	    if(!strcmp(selstage, "viability") || !strcmp(selstage, "none")){
	      pat_index = tmp; /* no fertility selection */
	      found = 1;
	    }
	    else if(gsl_rng_uniform(r) < gsl_vector_get(grid[source - 1].fitness, tmp) ){
	      pat_index = tmp;
	      found = 1;
	    }
	  }
	}

	mat_index = mat_index * kstrideInd;
	pat_index = pat_index * kstrideInd;


	/* progeny dispersal, identify (potential) sink deme */
	/* ------------------------------------------------- */
	sink = -1;

	if (!strcmp(migstage, "pollen")){ /* no progeny migration */
	  sink = k+1;
	}
	else { /* progeny migration */
	  gsl_ran_multinomial(r, ndemes+2, 1, migration->data, sinkvec->data);

	  /* sink will have zeros and one 1 */
	  for (m=0; m<ndemes+2; m++){
	    if (gsl_vector_uint_get(sinkvec,m) == 1){
	      sink = m;
	    }
	  }
	}



	/* make zygote, test zygote (relative fitness, viability) */
	/* ------------------------------------------------------ */
	for(homolog=0; homolog<kstrideInd; homolog = homolog + 2){ /* go through chromosome sets */

	  /* ------------ PATERNAL ---- parent ----------------- */
	  if(source == 0){ /* purebred immigrant; make new struct chromosome for it */
	    make_parental_gamete(&paternalgamete, 'a'); /* all alleles are 'a' */
	  }
	  else if(source == ndemes+1){ /* purebred immigrant; make new struct chromosome for it */
	    make_parental_gamete(&paternalgamete, 'b'); /* all alleles are 'b' */
	  }
	  else { /* go into correct deme; make gamete with recombi */
	    recombine(&grid[source - 1].sub_cur_generation[pat_index + homolog],
		      &grid[source - 1].sub_cur_generation[pat_index + homolog + 1],
		      &paternalgamete);
	  }

	  /* -------------- MATERNAL --------- parent --------------- */
	  /* moms are always residents; make gamete with recombi */
	  recombine(&grid[k].sub_cur_generation[mat_index + homolog],
		    &grid[k].sub_cur_generation[mat_index + homolog + 1],
		    &maternalgamete);

	  /* store gametes for each chromosome, first from mom, second from dad */
	  if(!strcmp(selstage, "both") || !strcmp(selstage, "viability")){
	    /* viability selection on potential progeny */
	    /* use storehomologViab instead of placezygote*/
	    init_chromosome(&potentialzygote[homolog], maternalgamete.nblocks);
	    storehomologViab(&potentialzygote[homolog], &maternalgamete);
	    init_chromosome(&potentialzygote[homolog+1], paternalgamete.nblocks);
	    storehomologViab(&potentialzygote[homolog+1], &paternalgamete);
	  }
	  else if (!strcmp(selstage, "fertility") || !strcmp(selstage, "none")){
	    /* store gamete in dyn_chromosome (nprocap is decremented below) */
	    if (sink == 0 || sink == ndemes+1){
	      ; /* dispersal into parental demes */
	    }
	    else if (sink > 0 && sink <= ndemes){
	      placezygote(maternalgamete, grid[sink - 1].dummyHead_sub_next_gen);
	      placezygote(paternalgamete, grid[sink - 1].dummyHead_sub_next_gen);
	    }
	    else {
	      fprintf(stderr, "Error in assigning sink deme for progeny\n");
	      exit(1);
	    }
	  }
	  else{
	    fprintf(stderr, "Please check spelling on selstage of selection:\"%s\"\n", selstage);
	    exit(1);
	  }
	  /* need to free memory temporarily allocated to subparts of gamete structs
	     in recombine and make_parental_gamete functions */

	  gsl_vector_char_free(maternalgamete.allele);
	  gsl_vector_float_free(maternalgamete.mapposition);
	  gsl_vector_char_free(paternalgamete.allele);
	  gsl_vector_float_free(paternalgamete.mapposition);

	} /* end loop through chromosomes */


	/* check fitness of zygote */
	if(!strcmp(selstage, "both") || !strcmp(selstage, "viability")){


	  /* compute zygote fitness */
	  if (nallselloci > 0){
	    getselgeno(geno, potentialzygote, 0,
		       epi, sel, nallselloci,
		       nepiloci, nselloci);
	    zygotefitness = getfitness(geno, epi,
				       nepiloci,
				       sel, nselloci, nallselloci,
				       selstrength,
				       grid[k].environment, grid[k].envselstrength);
	  }
	  else {
	    zygotefitness = 1.0;
	  }

	  if(gsl_rng_uniform(r) < zygotefitness){ /* viable zygote */
	    for(homolog=0; homolog<kstrideInd; homolog++){
	      /* copy zygote into next_gen */
	      /* store gamete in dyn_chromosome */
	      if (sink == 0 || sink == ndemes+1){
		; /* dispersal into parental demes */
	      }
	      else if (sink > 0 && sink <= ndemes){
		placezygote(potentialzygote[homolog], grid[sink - 1].dummyHead_sub_next_gen);
	      }
	      else {
		fprintf(stderr, "Error in assigning sink deme for progeny\n");
		exit(1);
	      }
	      /* free potentialzygote allocation for this homolog */
	      gsl_vector_char_free(potentialzygote[homolog].allele);
	      gsl_vector_float_free(potentialzygote[homolog].mapposition);
	    }
	    nprocap--;
	    if (sink > 0 && sink <= ndemes){
	      (grid[sink - 1].Npr)++;
	    }
	  }
	  else{
	    /* no viable zygote produced, don't decrement popcapacity */
	    for(homolog=0; homolog<kstrideInd; homolog++){
	      gsl_vector_char_free(potentialzygote[homolog].allele);
	      gsl_vector_float_free(potentialzygote[homolog].mapposition);
	    }
	  }
	}
	else{ /* no viability selection, thus progeny survives and decrements capacity */
	  nprocap--;
	  if (sink > 0 && sink <= ndemes){
	    (grid[sink - 1].Npr)++;
	  }
	} /* end homologs loop */

	/* nprocap--; only want to decrement if zygote viable */
      } /* end make zygotes loop */

      gsl_vector_uint_free(gameteID);

    }

    /* with progeny dispersal, demes might receive migrants from parental source */
    /* for this, the parental source sends as many migrants as a deme with similar
       characteristics as the peripheral deme, for convenience
       (i.e. potential #progeny given cap_ad, fecundity, and cap_pro (full fitness/fertility)) */

    if (!strcmp(migstage, "progeny") || !strcmp(migstage, "both")){


      /* immigrants from source 0 */
      p = gsl_vector_get(migration, 0);
      if (p > 0){ /* there is migration */
	/* identify number of migrants with pop characteristics from first deme */
	nproposs = 0;
	for (i=0; i<grid[0].cap_ad; i++){
	  tmpnpro = gsl_ran_poisson(r, MeanFecund);
	  nproposs += tmpnpro;
	}
	if (nproposs > grid[0].cap_pr) {
	  nproposs = grid[0].cap_pr;
	}
	allele = 'a';
	nproposs = gsl_ran_binomial(r, p, nproposs);

	/* add them to dyn_chromosome */
	make_parental_gamete(&maternalgamete, allele);

	for (i=0;i<nproposs;i++){
	  for (homolog=0; homolog<kstrideInd; homolog++){ /* go through chromosome sets */
	    placezygote(maternalgamete, grid[k].dummyHead_sub_next_gen);
	  }
	  (grid[k].Npr)++;

	}
	gsl_vector_char_free(maternalgamete.allele);
	gsl_vector_float_free(maternalgamete.mapposition);
      }


      /* immigrants from source 1 */
      p = gsl_vector_get(migration, ndemes+1);
      if (p > 0){ /* there is migration */
	/* identify number of migrants with pop characteristics from last deme */
	nproposs = 0;
	for (i=0; i<grid[ndemes-1].cap_ad; i++){
	  tmpnpro = gsl_ran_poisson(r, MeanFecund);
	  nproposs += tmpnpro;
	}
	if (nproposs > grid[ndemes-1].cap_pr) {
	  nproposs = grid[ndemes-1].cap_pr;
	}
	allele = 'b';
	nproposs = gsl_ran_binomial(r, p, nproposs);

	/* add them to dyn_chromosome */
	make_parental_gamete(&maternalgamete, allele);

	for (i=0;i<nproposs;i++){
	  for (homolog=0; homolog<kstrideInd; homolog++){ /* go through chromosome sets */
	    placezygote(maternalgamete, grid[k].dummyHead_sub_next_gen);
	  }
	  (grid[k].Npr)++;

	}
	gsl_vector_char_free(maternalgamete.allele);
	gsl_vector_float_free(maternalgamete.mapposition);
      }
    }

  } /* end deme loop */

  gsl_vector_uint_free(sourcevec);
  gsl_vector_uint_free(sinkvec);
  gsl_vector_free(migration_copy);
  gsl_vector_free(migration);
  if(nallselloci > 0){
    gsl_vector_int_free(geno);
  }

}


/* -------------------------- */
/* more simulation essentials */
/* -------------------------- */

void placezygote(struct chromosome gamete, struct dyn_chromosome * dummyHead){
  struct dyn_chromosome * newPtr;
  /* store gamete in dyn_chromosome */
  newPtr = malloc(sizeof(struct dyn_chromosome));
  if ( !newPtr ) {
    fprintf (stderr, "No memory available.");
    exit (-1);
  }

  init_dyn_chromosome(newPtr, gamete.nblocks); /* allocate memory */
  storehomolog(newPtr, &gamete); /* copy gamete into newPtr */
  /* insert newPtr into list; insertion is always just after the
       dummy/first pointer, end of list is indicated by a NULL nextPtr
       in the final field, this is a FILO list, which has to be
       reverse when copied to cur_generation */
  newPtr->nextPtr=dummyHead->nextPtr;  /* this will be NULL if newPtr is first item added */
  dummyHead->nextPtr = newPtr;
}



void copynext2cur(struct deme * grid){

  unsigned int k, missing, tmpN, n, chroms;
  long int m; /* long int is overkill but cannot use unsigned because of test m >= 0 below */
  gsl_vector_uint * proID;

  struct dyn_chromosome * nodePtr;


  for (k=0; k<ndemes; k++){

    /* to copy progeny to next gen at random, copy chromosomes as block for each ind */
    /* need to know who to copy */

    if (grid[k].Npr == 0){ /* nothing to copy */
      grid[k].N = 0;
    }
    else {
      proID = gsl_vector_uint_calloc(grid[k].Npr);
      if (grid[k].Npr > grid[k].cap_ad) { /* have to regulate to carrying capacity */

	gsl_vector_uint_set_all(proID, 0);
	for (m=0; m < grid[k].cap_ad; m++){
	  gsl_vector_uint_set(proID, m, 1);
	}
	gsl_ran_shuffle(r, proID->data, grid[k].Npr, sizeof(unsigned int));

	missing = 0;
      }
      else { /* don't have to regulate */
	gsl_vector_uint_set_all(proID, 1);
	missing = grid[k].cap_ad - grid[k].Npr;
      }

      missing *= kstrideInd;
      chroms = (grid[k].Npr * kstrideInd) - 1 ; /* need for switching to next ind */
      tmpN = 0;


      /* need to reverse order in the copy because next_gen has the last
	 entry at first spot in the list (after head), so counter descends */
      for(nodePtr = grid[k].dummyHead_sub_next_gen->nextPtr,
	    m=grid[k].gChromoNpop-1-missing; nodePtr != NULL && m >= 0;
	  nodePtr = nodePtr->nextPtr, chroms--){

	/* only copy if proID == 1 */
	if (gsl_vector_uint_get(proID, tmpN) == 1){

	  init_chromosome(&grid[k].sub_cur_generation[m], nodePtr->nblocks);

	  gsl_vector_char_memcpy(grid[k].sub_cur_generation[m].allele,
				 nodePtr->allele);
	  gsl_vector_float_memcpy(grid[k].sub_cur_generation[m].mapposition,
				 nodePtr->mapposition);
	  grid[k].sub_cur_generation[m].nblocks = nodePtr->nblocks;
	  m--;
	}
	else {
	  ; /* do nothing, but for-loop goes on */
	}
	if (chroms % kstrideInd == 0){
	  tmpN++;
	}

      }

      if(missing > 0){
        for(n=grid[k].gChromoNpop-missing;n < grid[k].gChromoNpop; n++){
      	/* linked list had fewer inds than carrying capacity of pop,
      	   fill the remained of the pop with zeros for nblocks */
      	grid[k].sub_cur_generation[n].nblocks = 0;
        }
      }


      grid[k].Npr = 0;
      grid[k].N = grid[k].cap_ad - missing/kstrideInd;

      gsl_vector_uint_free(proID);

    }
    gsl_vector_set_zero(grid[k].fitness);
    gsl_matrix_set_zero(grid[k].stats); /* clean out old stats */
    gsl_vector_set_zero(grid[k].junctions);
  }
}



double calc_heterozygosity(const struct chromosome * first,
			   const struct chromosome * second){
  /* calculate the fraction of the chromosome that is heterozygous  */
  int i, k, ctr, njunctions;
  double hetsum;
  char a1, a2;

  gsl_vector_float * alljunctions;
  gsl_vector_float_view all_unique_junctions; /* subvector */

  njunctions = first->nblocks + second->nblocks;
  /* find positions of all junctions for both homologs */
  alljunctions = gsl_vector_float_alloc(njunctions + 1);
  /* +1 spot for beginning of chromosome */

  /* concatenate the vectors, add in the range breaks, drop the upperbound */
  gsl_vector_float_set(alljunctions, 0, 0); /* set position 0 in array to zero */
  for(i=0; i<first->nblocks; i++){
    gsl_vector_float_set(alljunctions, i+1, gsl_vector_float_get(first->mapposition, i));
  }
  for(i=0; i<second->nblocks;i++){
    gsl_vector_float_set(alljunctions, i+1+first->nblocks,
			 gsl_vector_float_get(second->mapposition, i));
  }
  /*  gsl_vector_set(alljunctions, njunctions, 1); */ /* add on the end of the chromosome */
  gsl_sort_vector_float(alljunctions);

  /* remove redundant junctions */
  ctr=njunctions;
  for(k=0, i=1; i<ctr; i++){
    if(gsl_vector_float_get(alljunctions, i) == gsl_vector_float_get(alljunctions, i-1)){
      k++; /* vector offset to move items into correct position*/
      njunctions--;
    }
    gsl_vector_float_set(alljunctions, i-k, gsl_vector_float_get(alljunctions, i) );
  }
  /* cut end of vector */
  all_unique_junctions = gsl_vector_float_subvector(alljunctions, 0, njunctions);

  /* finally, truly calculate heterozygosity */
  for(i=0, k=0, ctr=0, hetsum=0;
      ctr<njunctions-1; /* -1 to avoid to step beyond chromosome end with ctr+1 */
      ctr++){
    a1 = gsl_vector_char_get(first->allele, i);
    a2 = gsl_vector_char_get(second->allele, k);
    if(tolower(a1) != tolower(a2)){ /* alleles can be 'a' or 'A', or 'b' or 'B' */
      hetsum = hetsum+ (gsl_vector_float_get(&all_unique_junctions.vector, ctr+1) -
			gsl_vector_float_get(&all_unique_junctions.vector, ctr));
      /* steps along unique junctions */
    }
    if( gsl_vector_float_get(first->mapposition, i) ==
	gsl_vector_float_get(&all_unique_junctions.vector, ctr+1) ){
      i++; /* move to next junction on first chromosome */
    }
    if( gsl_vector_float_get(second->mapposition, k) ==
	gsl_vector_float_get(&all_unique_junctions.vector, ctr+1) ){
      k++; /* move to next junction on second chromosome */
    }
  }

  gsl_vector_float_free(alljunctions);
  return(hetsum);
}



void copy_homolog (struct chromosome * gamete, const struct chromosome * parent){
  gamete->nblocks = parent->nblocks;
  gsl_vector_char_memcpy(gamete->allele, parent->allele);
  gsl_vector_float_memcpy(gamete->mapposition, parent->mapposition);
}




/* free memory allocated in cur_generation  */
void free_cur_gen(struct deme * grid){
  int k;
  unsigned int m;

  for (k=0; k<ndemes; k++){

    for (m=0; m < (grid[k].N * kstrideInd); m++){
      if(grid[k].sub_cur_generation[m].nblocks>0){
	grid[k].sub_cur_generation[m].nblocks = 0;
	gsl_vector_char_free(grid[k].sub_cur_generation[m].allele);
	gsl_vector_float_free(grid[k].sub_cur_generation[m].mapposition);
      }
    }
    /* some testing that no other chroms accidently have data */
    for (m = (grid[k].N * kstrideInd); m < grid[k].gChromoNpop; m++){
      if(grid[k].sub_cur_generation[m].nblocks>0){
	fprintf(stderr,"Found some unexpected data while freeing current generation\n");
	grid[k].sub_cur_generation[m].nblocks = 0;
	gsl_vector_char_free(grid[k].sub_cur_generation[m].allele);
	gsl_vector_float_free(grid[k].sub_cur_generation[m].mapposition);
      }
    }
  }
}



void dynamic_grid_free(struct deme * grid){

  int k;

  for (k=0; k<ndemes; k++){

    struct dyn_chromosome * nodePtr;
    struct dyn_chromosome * toremovePtr;

    for(nodePtr = grid[k].dummyHead_sub_next_gen->nextPtr; nodePtr != NULL; ){
      toremovePtr = nodePtr;
      nodePtr = toremovePtr->nextPtr;
      if(toremovePtr->nblocks > 0){
	gsl_vector_char_free(toremovePtr->allele);
	gsl_vector_float_free(toremovePtr->mapposition);
      }
      free(toremovePtr);
    }
  }
}

/* ------------------------------------------------------------------------- */
/* simple and typical version of summary, without storing size of all blocks */
void summarize_cur_gen(struct deme * grid, const int nallselloci,
		       const int nepiloci, const int nselloci, const double selstrength,
		       const struct epiloci * epi, const struct selloci * sel){
  int i, ind, k;
  double ablock, het, fit, njunct;
  gsl_vector_int * geno;

  if (nallselloci > 0){
    geno = gsl_vector_int_alloc(nallselloci);
  }

  /* for deme */
  for (k=0; k<ndemes; k++){

    if (grid[k].N > 0) {

      /* for individuals */
      for (i=0, ind=0; i < grid[k].N * kstrideInd; i = i+kstrideInd, ind++){
	/* i+kstrideInd goes to first chrom of next ind */
	summarize_ind(grid[k].sub_cur_generation, i, &ablock, &het, &njunct);
	gsl_matrix_set(grid[k].stats, ind, 0, ablock);
	gsl_matrix_set(grid[k].stats, ind, 1, het);
	gsl_vector_set(grid[k].junctions, ind, njunct);

	if (nallselloci > 0){
	  getselgeno(geno, grid[k].sub_cur_generation, ind,
		     epi, sel, nallselloci,
		     nepiloci, nselloci);

	  fit = getfitness(geno, epi,
			   nepiloci,
			   sel, nselloci, nallselloci,
			   selstrength,
			   grid[k].environment, grid[k].envselstrength);
	  gsl_vector_set(grid[k].fitness, ind, fit);

	}
	else{ /* no selected loci, full fitness */
	  gsl_vector_set(grid[k].fitness, ind, 1.0);
	}

      }
    }
  }

  if (nallselloci > 0){
    gsl_vector_int_free(geno);
  }

}



/*  summary for single ind  */
void summarize_ind(struct chromosome * ind, int offset, double * ablock, double * het,
		   double * njunct){
  int j,k, tmpnjunctions;
  char a1, a0;

  /* ablock is equivalent to a hybrid index.  It is the proportion of
     the genome from the 'a' species */

  *ablock = *het =  0;
  *njunct = 0.0;
  for(k=0; k<kstrideInd; k++){ /*  move along chromosomes of 1 ind */
    if(ind[offset+k].nblocks == 1 ){
      if(gsl_vector_char_get(ind[offset+k].allele, 0) == 'a' ||
	 gsl_vector_char_get(ind[offset+k].allele, 0) == 'A'){
	*ablock = *ablock + 1;
      }
    }
    else if (ind[offset+k].nblocks > 1){
      tmpnjunctions = ind[offset+k].nblocks - 1; /* don't count block that is not junction */
      for(j=0; j<ind[offset+k].nblocks; j++){ /* move along blocks */
	if(j==0){
	  /* take care of special case, the first block */
	  if(gsl_vector_char_get(ind[offset+k].allele, 0) == 'a' ||
	     gsl_vector_char_get(ind[offset+k].allele, 0) == 'A'){
	    /* proportion of chromosome with 'a' or 'A' ancestry in first block */
	    *ablock = *ablock + gsl_vector_float_get(ind[offset+k].mapposition, 0);
	  }
	}
	else{
	  if(gsl_vector_char_get(ind[offset+k].allele, j) == 'a' ||
	     gsl_vector_char_get(ind[offset+k].allele, j) == 'A'){
	    *ablock = *ablock +
	      gsl_vector_float_get(ind[offset+k].mapposition, j) -
	      /* size of next block with 'a' or 'A' ancestry */
	      gsl_vector_float_get(ind[offset+k].mapposition, j-1);
	  }
	  /* simplify number of junctions */
	  a1 = gsl_vector_char_get(ind[offset+k].allele, j);
	  a0 = gsl_vector_char_get(ind[offset+k].allele, j-1);
	  if(tolower(a1) == tolower(a0)) {
	    tmpnjunctions--;
	  }
	}
      }
      *njunct = *njunct + (double) tmpnjunctions;
    }
    if(0 == k%2){ /* for all even chromosomes */
      *het = *het + calc_heterozygosity(&ind[offset+k], &ind[offset+k+1]);
    }
  }
  *het = *het/khomologsN; /*  average fraction of the genome that is
			    heterozygous, across all chromosomes in an
			    individual */
  *ablock = *ablock / kstrideInd;

}


/* ------------------------------------ */
/* compute and write some summary stats */
/* ------------------------------------ */

/* one line per generation and deme, extra line for deme-unspecific stats */
void write_stats(FILE * file, int rep, int gen, struct deme * grid){

  int startdeme = 0, enddeme = ndemes;
  int k, m = 0, n = 0;
  double qall = 0.0, obshetall = 0.0, exphetall;
  double fst, fis, fisdeme, hetperdeme = 0.0, q, obshet, exphet;
  gsl_vector * qtmp;
  gsl_vector * hettmp;
  gsl_vector * juncttmp;

  /* fst is Nei's Fst = (Ht - Hs)/Ht */
  /* fis is (He - Ho)/He */

  if(dispersalmode==1){ /* omit first and last deme */
    startdeme=1;
    enddeme=ndemes-1;
  }

  for (k=startdeme; k<enddeme; k++){

    if (grid[k].N > 0) {

      /* within deme stats */
      qtmp = gsl_vector_alloc(grid[k].cap_ad);
      hettmp = gsl_vector_alloc(grid[k].cap_ad);
      juncttmp = gsl_vector_alloc(grid[k].cap_ad);

      gsl_matrix_get_col(qtmp, grid[k].stats, 0);
      gsl_matrix_get_col(hettmp, grid[k].stats, 1);
      gsl_vector_memcpy(juncttmp, grid[k].junctions);

      q = gsl_stats_mean(qtmp->data, 1, grid[k].N);
      obshet = gsl_stats_mean(hettmp->data, 1, grid[k].N);
      exphet = 2 * q * (1-q);
      fisdeme = (exphet - obshet)/exphet;

      /* metapopulation stats */
      qall += q * grid[k].N;
      obshetall += obshet * grid[k].N;
      hetperdeme += exphet;
      n += grid[k].N; /* ind counter */
      m++; /* deme counter */


      fprintf(file, "%d,%d,%d,%d,%f,%f,%f,%f,%f,%f,%f,%f,%f,NA,NA\n",
	      rep+1,gen,(dispersalmode==1)?k:k+1,grid[k].N,q,gsl_stats_sd(qtmp->data, 1, grid[k].N),
	      obshet,gsl_stats_sd(hettmp->data, 1, grid[k].N),
	      gsl_stats_mean((grid[k].fitness)->data, 1, grid[k].N),
	      gsl_stats_sd((grid[k].fitness)->data, 1, grid[k].N),
	      gsl_stats_mean(juncttmp->data, 1, grid[k].N),
	      gsl_stats_sd(juncttmp->data, 1, grid[k].N),
	      fisdeme);


      gsl_vector_free(qtmp);
      gsl_vector_free(hettmp);
      gsl_vector_free(juncttmp);

    }
    else { /* deme empty */
      fprintf(file, "%d,%d,%d,%d,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA\n",
	      rep+1,gen,(dispersalmode==1)?k:k+1,grid[k].N);
    }

  }

  /* metapopulation stats */
  /* (for dispersalmode==1 these things should have been calculated for middle demes only above) */
  qall = qall/n;
  obshetall = obshetall/n;
  exphetall = 2 * qall * (1 - qall);
  hetperdeme = hetperdeme/m;

  fis = (exphetall - obshetall)/exphetall;
  fst = (exphetall - hetperdeme)/exphetall;

  fprintf(file, "%d,%d,meta,%d,%f,NA,%f,NA,NA,NA,NA,NA,NA,%f,%f\n",
	  rep+1,gen,n,qall,
	  obshetall,fis,fst);

}


/* compute and write some summary stats for selected loci */
/* one line per replicate, generation and deme */
void write_selstats(FILE * file, const int rep, const int gen, struct deme * grid,
		    const struct epiloci * epi, const struct selloci * sel,
		    const int nepiloci, const int nselloci, const int nallselloci) {

  int startdeme = 0, enddeme = ndemes;
  int ind, k, geno, j, l;
  double het, af;

  if(dispersalmode==1){ /* omit first and last deme */
    startdeme=1;
    enddeme=ndemes-1;
  }

  /* for deme */
  for (k=startdeme; k<enddeme; k++){

    /* print replicate, generation, deme */
    fprintf(file, "%d,%d,%d", rep+1, gen, (dispersalmode==1)?k:k+1);

    /* for selected loci */
    for(j=0, l=0; j<nallselloci; j++){
      het = 0.0;
      af = 0.0;
      /* go through epiloci first, then through other loci */
      if (j<nepiloci){
	if (grid[k].N > 0) { /* there are data */
	  for(ind=0; ind < grid[k].N; ind++){
	    geno = getgenotype(grid[k].sub_cur_generation, ind,
			       gsl_vector_int_get(epi->chr, j),
			       gsl_vector_get(epi->posdecimal, j));
	    af += (double) geno;
	    het += ((geno==1) ? 1.0 : 0.0);
	  }
	  af = af/(2*grid[k].N);
	  het = het/grid[k].N;
	  fprintf(file, ",%f,%f",af,het);
	}
	else { /* deme empty */
	  fprintf(file, ",NA,NA");
	}
      }
      else if (l<nselloci){
	if (grid[k].N > 0) { /* there are data */
	  for(ind=0; ind < grid[k].N; ind++){
	    geno = getgenotype(grid[k].sub_cur_generation, ind,
			       gsl_vector_int_get(sel->chr, l),
			       gsl_vector_get(sel->posdecimal, l));
	    af += (double) geno;
	    het += ((geno==1) ? 1.0 : 0.0);
	  }
	  af = af/(2*grid[k].N);
	  het = het/grid[k].N;
	  fprintf(file, ",%f,%f",af,het);
	}
	else { /* deme empty */
	  fprintf(file, ",NA,NA");
	}

	l++;
      }
    }
    fprintf(file, "\n");

  }
}



int getgenotype(struct chromosome * cur_generation, int ind, int chr,
		 double position){
  int index, genotype;

  genotype = 0;
  /* grab allele data at given index from homolog */
  /* find block that contains selected locus */
  index = find_position(position, &cur_generation[(ind * kstrideInd) + chr * 2]) ;
  if(gsl_vector_char_get(cur_generation[(ind * kstrideInd) + chr * 2].allele, index) == 'a' ||
     gsl_vector_char_get(cur_generation[(ind * kstrideInd) + chr * 2].allele, index) == 'A'){
    genotype++;
  }
  index = find_position(position, &cur_generation[(ind * kstrideInd) + chr * 2 + 1]) ;
  if(gsl_vector_char_get(cur_generation[(ind * kstrideInd) + chr * 2 + 1].allele, index) == 'a' ||
     gsl_vector_char_get(cur_generation[(ind * kstrideInd) + chr * 2 + 1].allele, index) == 'A' ){
    genotype++;
  }
  return(genotype); /* returns number of 'a' or 'A' alleles at selected site */
}

/* get all selected sites and store in vector */
void getselgeno(gsl_vector_int * geno, struct chromosome * cur_generation, const int ind,
		const struct epiloci * epi, const struct selloci * sel, const int nallselloci,
		const int nepiloci, const int nselloci){

  int j, k;

  /* go through epiloci first, then through other loci */
  for(j=0, k=0; j<nallselloci; j++){
    if (j<nepiloci){
      gsl_vector_int_set(geno, j, getgenotype(cur_generation, ind,
					      gsl_vector_int_get(epi->chr, j),
					      gsl_vector_get(epi->posdecimal, j)));
    }
    else if (k<nselloci){
      gsl_vector_int_set(geno, j, getgenotype(cur_generation, ind,
					      gsl_vector_int_get(sel->chr, k),
					      gsl_vector_get(sel->posdecimal, k)));
      k++;
    }
  }

}


int gethaplotype(struct chromosome * cur_generation, int ind, int chr, int copy,
		 double position){
  int index, haplotype;

  haplotype = 0;
  /* grab allele data at given index from homolog */
  /* find block index that contains selected locus */
  index = find_position(position, &cur_generation[(ind * kstrideInd) + chr * 2 + copy]) ;
  if(gsl_vector_char_get(cur_generation[(ind * kstrideInd) + chr * 2 + copy].allele,
			 index) == 'a' ||
     gsl_vector_char_get(cur_generation[(ind * kstrideInd) + chr * 2 + copy].allele,
			 index) == 'A' ){
    haplotype++;
  }
  return(haplotype); /* returns number of 'a' or 'A' alleles at selected site */
}


/* same as gethaplotype, but return char instead */
char gethaplotype_char(struct chromosome * cur_generation, int ind, int chr, int copy,
		 double position){
  int index;
  char haplotype;

  /* grab allele data at given index from homolog */
  /* find block index that contains selected locus */
  index = find_position(position, &cur_generation[(ind * kstrideInd) + chr * 2 + copy]) ;
  haplotype = gsl_vector_char_get(cur_generation[(ind * kstrideInd) + chr * 2 + copy].allele,
				index);
  return(haplotype); /* returns 'a', 'A', 'b' or 'B' at locus */
}


double getdyefreq(struct chromosome * cur_generation, const int N, const int chr,
		const double position, const char dye){

  int ind, index;
  int dyecount = 0;

  for (ind = 0; ind < N; ind++){
    index = find_position(position, &cur_generation[(ind * kstrideInd) + chr * 2]) ;
    if(gsl_vector_char_get(cur_generation[(ind * kstrideInd) + chr * 2].allele,
			   index) == dye){
      dyecount++;
    }
    index = find_position(position, &cur_generation[(ind * kstrideInd) + chr * 2 + 1]) ;
    if(gsl_vector_char_get(cur_generation[(ind * kstrideInd) + chr * 2 + 1].allele,
			   index) == dye){
      dyecount++;
    }
  }

  return((double) dyecount / (double) (N * 2));
}


double getfitness(gsl_vector_int * genotype, const struct epiloci * epi,
		  const int nepiloci,
		  const struct selloci * sel, const int nselloci, const int nallselloci,
		  const double selstrength,
		  const gsl_vector * envi, const gsl_vector * envsel){

  int i, j=0, k;
  int g1, g2, g;
  double fit = 1.0, mismatch = 0.0;


  /* epiloci */
  /* ------- */

  /* in genotype, number of 'a' or 'A' alleles are stored;
     need to translate genotype in number of derived alleles */

  if(nepiloci > 0){

    /* for complex */
    for (i=0; i < epi->ncomplexes; i++){
      /* for nepiloci within complex */

      /* dmi */
      if (gsl_vector_char_get(epi->type, (i * epi->xplets)) == 'd'){
	for (j=0; j < (epi->xplets)-1; j++){
	  if (gsl_vector_char_get(epi->ancallele, (i * epi->xplets) + j) == 'a'){
	    g1 = 2 - gsl_vector_int_get(genotype, (i * epi->xplets) + j); /* # alleles not 'a' */
	  }
	  else if (gsl_vector_char_get(epi->ancallele, (i * epi->xplets) + j) == 'b'){
	    g1 = gsl_vector_int_get(genotype, (i * epi->xplets) + j); /* # alleles 'a' */
	  }
	  else {
	    fprintf(stderr, "Error: ancestral allele unknown\n");
	    exit(1);
	  }

	  /* get next locus genotype and calculate s */
	  for (k=j+1; k < epi->xplets; k++){
	    if (gsl_vector_char_get(epi->ancallele, (i * epi->xplets) + j) ==
		gsl_vector_char_get(epi->ancallele, (i * epi->xplets) + k)) {
	      ; /* no derived-derived incompatibility possible, fit remains unchanged */
	    }
	    else {
	      if (gsl_vector_char_get(epi->ancallele, (i * epi->xplets) + k) == 'a'){
		g2 = 2 - gsl_vector_int_get(genotype, (i * epi->xplets) + k); /* # alleles not 'a' */
	      }
	      else if (gsl_vector_char_get(epi->ancallele, (i * epi->xplets) + k) == 'b'){
		g2 = gsl_vector_int_get(genotype, (i * epi->xplets) + k); /* # alleles 'a' */
	      }
	      else {
		fprintf(stderr, "Error: ancestral allele unknown\n");
		exit(1);
	      }
	      fit = fit * (epi->dmimat)[g1][g2]; /* multiplicative */
	    }
	  }
	}
      }
      /* linear pathway */
      else if (gsl_vector_char_get(epi->type, (i * epi->xplets)) == 'p'){
	for (j=0; j < (epi->xplets)-1; j++){
	  g1 = gsl_vector_int_get(genotype, (i * epi->xplets) + j); /* # alleles 'a' */
	  g2 = gsl_vector_int_get(genotype, (i * epi->xplets) + j+1); /* # alleles 'a' */
	  fit = fit * (epi->pathmat)[g1][g2]; /* multiplicative */
	}
      }
      else {
	fprintf(stderr, "Error: selection type in epifile unknown\n");
	exit(1);
      }
    }

  }
  /* selloci */
  /* ------- */

  if (nselloci > 0){

    for (i=0, j=i+nepiloci; i<nselloci; i++, j++){
      /* i goes through selloci, j goes through genotypes */

      /* environmental selection */
      if (gsl_vector_char_get(sel->type, i) == 'e'){
	/* calculate mismatch between genotype and environment */
	/* multiply by selection factor */
	/* fitness is calculated multiplicatively here */
	g = gsl_vector_int_get(genotype, j);
	if(g == 2){/* 'a' homozygote */
	  mismatch = fabs(gsl_vector_get(envi, i) - gsl_vector_get(sel->valaa, i));
	}
	else if (g == 1){ /* heterozygote */
	  mismatch = fabs(gsl_vector_get(envi, i) - gsl_vector_get(sel->valab, i));
	}
	else if(g == 0){
	  mismatch = fabs(gsl_vector_get(envi, i) - gsl_vector_get(sel->valbb, i));
	}
	fit = fit * (1 - (mismatch * gsl_vector_get(envsel, i)));

      }
      /* overdominance */
      else if (gsl_vector_char_get(sel->type, i) == 'o'){
	g = gsl_vector_int_get(genotype, j);
	if(g == 2){/* 'a' homozygote */
	  fit = fit * (1 - selstrength);
	}
	else if (g == 1){ /* heterozygote */
	  ;
	}
	else if(g == 0){
	  fit = fit * (1 - selstrength);
	}
      }
      /* underdominance */
      else if (gsl_vector_char_get(sel->type, i) == 'u'){
	g = gsl_vector_int_get(genotype, j);
	if(g == 2){/* 'a' homozygote */
	  ;
	}
	else if (g == 1){ /* heterozygote */
	  fit = fit * (1 - selstrength);
	}
	else if(g == 0){
	  ;
	}
      }
      else {
	fprintf(stderr, "Error: selection type in selfile unknown\n");
	exit(1);
      }
    }

  }

  if (nepiloci+nselloci<nallselloci){
    fprintf(stderr, "Warning: problem with getting the number of selected loci right\n");
  }

  return(fit);
}





int make_parental_gamete(struct chromosome * gam, const char state){

  char altstate = 'c';
  if(stain == 1 && state == 'a'){
    altstate = 'A';
  }
  else if(stain == 1 && state == 'b'){
    altstate = 'B';
  }
  else {
    altstate = state;
  }


  init_chromosome(gam, 1);
  gsl_vector_char_set(gam->allele, 0, altstate);
  gsl_vector_float_set(gam->mapposition, 0, 1);
  gam->nblocks = 1;
  return(0);
}


/* ------ */
/* output */
/* ------ */

void print_file(FILE * file,
		int rep, int gen, gsl_vector * markerloc, int markers_per_chr,
		struct deme * grid){

  int startdeme = 0, enddeme = ndemes;
  int chr, g, geno, ind, k;

  if(dispersalmode==1){ /* omit first and last deme */
    startdeme=1;
    enddeme=ndemes-1;
  }

  for (k=startdeme; k<enddeme; k++){

    for(ind=0; ind < grid[k].cap_ad; ind++){
      fprintf(file, "%d,%d,%d,%d",
	      rep+1, gen, (dispersalmode==1)?k:k+1, ind+1 ); /* replicate, generation, deme, individual */

      if (ind < grid[k].N){ /* there are data (deme might not be filled to capacity) */
	/* hybridindex, heterozygosity, fitness, # junctions */
	fprintf(file, ",%f", gsl_matrix_get(grid[k].stats,ind,0));
	fprintf(file, ",%f", gsl_matrix_get(grid[k].stats,ind,1));
	fprintf(file, ",%f", gsl_vector_get(grid[k].fitness,ind));
	fprintf(file, ",%.0f", gsl_vector_get(grid[k].junctions,ind));
	/* genotype */
	for(chr=0; chr < khomologsN; chr++){
	  for(g=0; g<markers_per_chr; g++){
	    geno = getgenotype(grid[k].sub_cur_generation, ind, chr,
			       gsl_vector_get(markerloc, g));
	    fprintf(file, ",%d", geno);
	  }
	}
	/* get line ending correct */
	fprintf(file, "\n");
      }

      else { /* not filled to carrying capacity */
	fprintf(file, ",NA,NA,NA,NA"); /* no hybridindex, etc. */
	for(g=0; g < markers_per_chr * khomologsN; g++){
	  fprintf(file, ",NA");
	}
	fprintf(file, "\n");
      }
    }
  }
}


void print_haplofile(FILE * file,
		int rep, int gen, gsl_vector * markerloc, int markers_per_chr,
		     struct deme * grid){

  int startdeme = 0, enddeme = ndemes;
  int chr, g, haplo, ind, k, copy;

  if(dispersalmode==1){ /* omit first and last deme */
    startdeme=1;
    enddeme=ndemes-1;
  }

  for (k=startdeme; k<enddeme; k++){

    for(ind=0; ind < grid[k].cap_ad; ind++){
      for(copy=0; copy<2; copy++){

	fprintf(file, "%d,%d,%d,%d",
		rep+1, gen, (dispersalmode==1)?k:k+1, ind+1 ); /* replicate, generation, deme, individual */

	if (ind < grid[k].N){ /* there are data (deme might not be filled to capacity) */
	  /* hybridindex, heterozygosity, fitness, # junctions  */
	  fprintf(file, ",%f", gsl_matrix_get(grid[k].stats,ind,0));
	  fprintf(file, ",%f", gsl_matrix_get(grid[k].stats,ind,1));
	  fprintf(file, ",%f", gsl_vector_get(grid[k].fitness,ind));
	  fprintf(file, ",%.0f", gsl_vector_get(grid[k].junctions,ind));
	  /* genotype */
	  for(chr=0; chr < khomologsN; chr++){
	    for(g=0; g<markers_per_chr; g++){
	      haplo = gethaplotype(grid[k].sub_cur_generation, ind, chr, copy,
				 gsl_vector_get(markerloc, g));
	      fprintf(file, ",%d", haplo);
	    }
	  }
	  /* get line ending correct */
	  fprintf(file, "\n");
	}

	else { /* not filled to carrying capacity */
	  fprintf(file, ",NA,NA,NA,NA"); /* no hybridindex, etc. */
	  for(g=0; g < markers_per_chr * khomologsN; g++){
	    fprintf(file, ",NA");
	  }
	  fprintf(file, "\n");
	}
      }
    }
  }
}


/* same as print_haplofile, but print alleles as chars instead */
void print_haplofile_char(FILE * file,
		int rep, int gen, gsl_vector * markerloc, int markers_per_chr,
			  struct deme * grid){

  int startdeme = 0, enddeme = ndemes;
  int chr, g, ind, k, copy;
  char haplo;

  if(dispersalmode==1){ /* omit first and last deme */
    startdeme=1;
    enddeme=ndemes-1;
  }

  for (k=startdeme; k<enddeme; k++){

    for(ind=0; ind < grid[k].cap_ad; ind++){
      for(copy=0; copy<2; copy++){

	fprintf(file, "%d,%d,%d,%d",
		rep+1, gen, (dispersalmode==1)?k:k+1, ind+1 ); /* replicate, generation, deme, individual */

	if (ind < grid[k].N){ /* there are data (deme might not be filled to capacity) */
	  /* hybridindex, heterozygosity, fitness, # junctions  */
	  fprintf(file, ",%f", gsl_matrix_get(grid[k].stats,ind,0));
	  fprintf(file, ",%f", gsl_matrix_get(grid[k].stats,ind,1));
	  fprintf(file, ",%f", gsl_vector_get(grid[k].fitness,ind));
	  fprintf(file, ",%.0f", gsl_vector_get(grid[k].junctions,ind));
	  /* genotype */
	  for(chr=0; chr < khomologsN; chr++){
	    for(g=0; g<markers_per_chr; g++){
	      haplo = gethaplotype_char(grid[k].sub_cur_generation, ind, chr, copy,
				 gsl_vector_get(markerloc, g));
	      fprintf(file, ",%c", haplo);
	    }
	  }
	  /* get line ending correct */
	  fprintf(file, "\n");
	}

	else { /* not filled to carrying capacity */
	  fprintf(file, ",NA,NA,NA,NA"); /* no hybridindex, etc. */
	  for(g=0; g < markers_per_chr * khomologsN; g++){
	    fprintf(file, ",NA");
	  }
	  fprintf(file, "\n");
	}

      }
    }
  }
}



void print_stainfile(FILE * file,
		     int rep, int gen, gsl_vector * markerloc, int markers_per_chr,
		     struct deme * grid){

  int chr, g, i, k;
  double dyefreq;
  char dye[] = {'A','B'};

  for (k=0; k<ndemes; k++){

    for(i=0;i<2;i++){ /* go through dyes */
      fprintf(file, "%d,%d,%d,%c",
	      rep+1, gen, k+1, dye[i]); /* replicate, generation, deme, dye */

      if (grid[k].N > 0){ /* there are data (deme might not be filled to capacity) */
	/* get dyefreq */
	for(chr=0; chr < khomologsN; chr++){
	  for(g=0; g<markers_per_chr; g++){
	    dyefreq = getdyefreq(grid[k].sub_cur_generation, grid[k].N, chr,
				 gsl_vector_get(markerloc, g), dye[i]);
	    fprintf(file, ",%f", dyefreq);
	  }
	}
	/* get line ending correct */
	fprintf(file, "\n");
      }

      else { /* not filled to carrying capacity */
	for(g=0; g < markers_per_chr * khomologsN; g++){
	  fprintf(file, ",NA");
	}
	fprintf(file, "\n");
      }
    }
  }
}



void print_selfilehead(FILE * file,
		       const struct epiloci * epi, const struct selloci * sel,
		       const int nepiloci, const int nselloci){

  int i;

  /* print info about selected loci */
  fprintf(file, "locuscount");
  for (i = 0; i < nepiloci; i++){
    fprintf(file, ",%d", gsl_vector_int_get(epi->locuscount, i));
  }
  for (i = 0; i < nselloci; i++){
    fprintf(file, ",%d", gsl_vector_int_get(sel->locuscount, i));
  }
  fprintf(file, "\ncomplex");
  for (i = 0; i < nepiloci; i++){
    fprintf(file, ",%d", gsl_vector_int_get(epi->complex, i));
  }
  for (i = 0; i < nselloci; i++){
    fprintf(file, ",%d", -1);
  }
  fprintf(file, "\ncompcount");
  for (i = 0; i < nepiloci; i++){
    fprintf(file, ",%d", gsl_vector_int_get(epi->compcount, i));
  }
  for (i = 0; i < nselloci; i++){
    fprintf(file, ",%d", -1);
  }
  fprintf(file, "\ntype");
  for (i = 0; i < nepiloci; i++){
    fprintf(file, ",%c", gsl_vector_char_get(epi->type, i));
  }
  for (i = 0; i < nselloci; i++){
    fprintf(file, ",%c", gsl_vector_char_get(sel->type, i));
  }
  fprintf(file, "\nchr");
  for (i = 0; i < nepiloci; i++){
    fprintf(file, ",%d", gsl_vector_int_get(epi->chr, i));
  }
  for (i = 0; i < nselloci; i++){
    fprintf(file, ",%d", gsl_vector_int_get(sel->chr, i));
  }
  fprintf(file, "\nposdecimal");
  for (i = 0; i < nepiloci; i++){
    fprintf(file, ",%.4f", gsl_vector_get(epi->posdecimal, i));
  }
  for (i = 0; i < nselloci; i++){
    fprintf(file, ",%.4f", gsl_vector_get(sel->posdecimal, i));
  }
  fprintf(file, "\nancallele");
  for (i = 0; i < nepiloci; i++){
    fprintf(file, ",%c", gsl_vector_char_get(epi->ancallele, i));
  }
  for (i = 0; i < nselloci; i++){
    fprintf(file, ",%c", 'w');
  }
  fprintf(file, "\nvalaa");
  for (i = 0; i < nepiloci; i++){
    fprintf(file, ",%.0f", -1.0);
  }
  for (i = 0; i < nselloci; i++){
    fprintf(file, ",%.2f", gsl_vector_get(sel->valaa, i));
  }
  fprintf(file, "\nvalab");
  for (i = 0; i < nepiloci; i++){
    fprintf(file, ",%.0f", -1.0);
  }
  for (i = 0; i < nselloci; i++){
    fprintf(file, ",%.2f", gsl_vector_get(sel->valab, i));
  }
  fprintf(file, "\nvalbb");
  for (i = 0; i < nepiloci; i++){
    fprintf(file, ",%.0f", -1.0);
  }
  for (i = 0; i < nselloci; i++){
    fprintf(file, ",%.2f", gsl_vector_get(sel->valbb, i));
  }
  fprintf(file, "\n");

}


/* print selected haplotypes */
void print_selfile(FILE * file, struct deme * grid, const int rep, const int gen,
		   const struct epiloci * epi, const struct selloci * sel,
		   const int nepiloci, const int nselloci, const int nallselloci){

  int startdeme = 0, enddeme = ndemes;
  int ind, k, copy, j, l, haplo;

  if(dispersalmode==1){ /* omit first and last deme */
    startdeme=1;
    enddeme=ndemes-1;
  }

  /* for deme */
  for (k=startdeme; k<enddeme; k++){

    for(ind=0; ind < grid[k].cap_ad; ind++){
      for(copy=0; copy<2; copy++){

	/* print replicate, generation, deme, individual */
	fprintf(file, "r%d_g%d_d%d_i%d", rep+1, gen, (dispersalmode==1)?k:k+1, ind+1);

	if (ind < grid[k].N){ /* there are data (deme might not be filled to capacity) */

	  /* go through epiloci first, then through other loci */
	  for(j=0, l=0; j<nallselloci; j++){
	    if (j<nepiloci){
	      haplo = gethaplotype(grid[k].sub_cur_generation, ind,
				   gsl_vector_int_get(epi->chr, j), copy,
				   gsl_vector_get(epi->posdecimal, j));
	    }
	    else if (l<nselloci){
	      haplo = gethaplotype(grid[k].sub_cur_generation, ind,
				   gsl_vector_int_get(sel->chr, l), copy,
				   gsl_vector_get(sel->posdecimal, l));
	      l++;
	    }
	    fprintf(file, ",%d", haplo);
	  }

	}
	else { /* not filled to carrying capacity */
	  for (j=0; j<nallselloci; j++){
	    fprintf(file, ",NA");
	  }
	}
	fprintf(file, "\n");
      }
    }
  }

}


void license(char * name){
  fprintf(stderr, "\n\n%s, a program to simulate hybridization and admixture.\n", name);
  fprintf(stderr, "Copyright (C) 2015, Doro Lindtke and Alex Buerkle\n\n");
  fprintf(stderr, "This program is free software: you can redistribute it and/or modify\n");
  fprintf(stderr, "it under the terms of the GNU General Public License as published\n");
  fprintf(stderr, "by the Free Software Foundation, either version 3 of the License,\n");
  fprintf(stderr, "or any later version.\n\n");
  fprintf(stderr, "This program is distributed in the hope that it will be useful,\n");
  fprintf(stderr, "but WITHOUT ANY WARRANTY; without even the implied warranty of\n");
  fprintf(stderr, "MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the\n");
  fprintf(stderr, "GNU General Public License for more details.\n\n");
  fprintf(stderr, "You should have received a copy of the GNU General Public License\n");
  fprintf(stderr, "along with this program. If not, see <http://www.gnu.org/licenses/>.\n\n\n");
  exit(1);
}

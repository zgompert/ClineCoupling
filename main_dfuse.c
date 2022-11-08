/* File: main_dfuse.c */

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

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <getopt.h>
#include <string.h>
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

/* ##################################################################### */
/* -----compile:----- */
/* gcc -O2 -Wall -o dfuse main_dfuse.c func_dfuse.c -lgsl -lgslcblas -lm */
/* ##################################################################### */


/* variables declared in head_dfuse.h as extern */
gsl_rng * r;  /* global state variable for random number generator */
int ndemes; /* number of demes */
int stain; /* set to true at generation x when staining starts */
int dispersalmode; /*  0=infinite source (default) or 1=finite source */

int main (int argc, char *argv[]){

  gsl_rng_env_setup();
  r = gsl_rng_alloc (gsl_rng_default);
  srand(time(NULL));
  int rng_seed;

  FILE *fp; /* for output files */
  int rep, repN = 1, gen, finalgen = 10;
  int outgen = 10; /* print every x generation */
  int selstats = 0; /* stats for selected loci every generation (0 = no) */
  int staingen = 0; /* generation when staining starts (0 = no staining) */
  int markers_per_chr = 51;
  double selstrength = 0.0;
  double ridge = 0.0; /* selection against ancestral alleles along ridge */
  int recessivity = 0; /* recessivity for DMIs (0 = dominant) */
  int nallselloci; /* total number of selected loci */
  int nepiloci; /* number of epistatically selected loci */
  int nselloci, nenviloci; /* selected, non-epistatic loci;
			      environmentally selected loci */

  char selstage[25] = "viability"; /* viability, fertility, both, or none */
  char migstage[25] = "progeny"; /* pollen, progeny or both */
  char multadd[25] = "mult"; /* multiplicative or additive */
  /* >>>> additive fitness accumulation not implemented yet */

  char mainfile[64] = "testing"; /* base name; for main output file */
  char haplofile[64]; /* for haplotypes output file */
  char seloutfile[64]; /* for selected loci output file */
  char selstatsfile[64]; /* for selected loci statistics output file */
  char statsfile[64]; /* for statistics output file */
  char stainfile[64]; /* for output file of freqs of stained alleles */
  char logfile[64]; /* for settings info file */

  char mydemefile[100] = "undefined";
  char myepifile[100] = "undefined";
  char myselfile[100] = "undefined";
  char myenvifile[100] = "undefined";

  int ch = 0;
  int i, j, k;

  double mig = 0.05; /* migration rate */

  struct auxcont aux;
  struct deme * grid;
  struct epiloci epi;
  struct selloci sel;

  gsl_vector * markerloc;
  dispersalmode = 0; /* set default */

  /* ------------------ */
  /* optional arguments */
  /* ------------------ */
  /* d: demefile, e: epifile, s: selfile, v: envifile, o: outfile basename */
  /* r: replicates, g: generations, G: steps to print, O: selected loci stats */
  /* m: migration, c: sel coeff, a: selection against ancestral alleles */
  /* R: recessivity for DMIs, S: sel stage, M: migr stage, D: dispersal mode */
  /* l: marker loci, f: fitness accumul model, Y: generation to start staining */

  while ((ch = getopt(argc, argv, "d:e:s:v:o:r:g:G:O:m:c:a:R:S:M:D:l:f:Y:Lh")) != -1){
    switch(ch){
    case 'd':
      strcpy(mydemefile, optarg);
      break;
    case 'e':
      strcpy(myepifile, optarg);
      break;
    case 's':
      strcpy(myselfile, optarg);
      break;
    case 'v':
      strcpy(myenvifile, optarg);
      break;
    case 'o':
      strcpy(mainfile, optarg);
      break;
    case 'r':
      repN = atoi(optarg);
      break;
    case 'g':
      finalgen = atoi(optarg);
      break;
    case 'G':
      outgen = atoi(optarg);
      break;
    case 'O':
      selstats = atoi(optarg);
      break;
    case 'm':
      mig = atof(optarg);
      break;
    case 'c':
      selstrength = atof(optarg);
      break;
    case 'a':
      ridge = atof(optarg);
      break;
    case 'R':
      recessivity = atoi(optarg);
      break;
    case 'S':
      strcpy(selstage, optarg);
      break;
    case 'M':
      strcpy(migstage, optarg);
      break;
    case 'D':
      dispersalmode = atoi(optarg);
      break;
    case 'l':
      markers_per_chr = atoi(optarg);
      break;
    case 'f':
      strcpy(multadd, optarg);
      break;
    case 'Y':
      staingen = atoi(optarg);
      break;
    case 'L':
      license(argv[0]);
      break;
    case 'h':
    default:
      usage(argv[0]);
    }
  }

  /* running and license statement */
  fprintf(stderr,"\nStarting dfuse. Type dfuse -h to see options.\n");
  fprintf(stderr,"Copyright (C) 2015, Doro Lindtke and Alex Buerkle. Type dfuse -L for details.\n\n");


  /* output files */
  strcat(strcpy(haplofile,mainfile),".haplo");
  strcat(strcpy(seloutfile,mainfile),".sel");
  strcat(strcpy(selstatsfile,mainfile),".selstats");
  strcat(strcpy(statsfile,mainfile),".stats");
  strcat(strcpy(stainfile,mainfile),".stain");
  strcat(strcpy(logfile,mainfile),".log");
  strcat(mainfile,".main");


  markerloc = gsl_vector_calloc(markers_per_chr);
  for(i=0; i<markers_per_chr; i++){
    gsl_vector_set(markerloc, i, (double) i/(markers_per_chr-1));
  }


  /* ----------------------- */
  /* some tests and warnings */
  /* ----------------------- */
  if(!strcmp(selstage,"none") && (strcmp(myepifile,"undefined") ||
				  strcmp(myselfile,"undefined") ||
				  strcmp(myenvifile,"undefined") ||
				  selstrength > 0)) {
    fprintf(stderr,"Warning: selection stage is \"none\" but selection params provided\n");
  }
  if(strcmp(myepifile,"undefined") && selstrength == 0) {
    fprintf(stderr,"Warning: epistasis file provided but selection coeff is zero\n");
  }
  if(!strcmp(myepifile,"undefined") && ridge > 0) {
    fprintf(stderr,"Warning: selection along DMI ridge requires epistasis file\n");
  }
  if(strcmp(myselfile,"undefined") && selstrength == 0) {
    fprintf(stderr,"Warning: single-locus selection file provided but selection coeff is zero, \n");
    fprintf(stderr,"         verify that environmental selection coeff is specified in envifile\n");
  }
  if(!(!strcmp(selstage,"viability") || !strcmp(selstage,"fertility") ||
       !strcmp(selstage,"both") || !strcmp(selstage,"none"))) {
    fprintf(stderr, "Error: \"%s\" is no valid stage of selection\n", selstage);
    exit(1);
  }
  if(!(!strcmp(migstage,"pollen") || !strcmp(migstage,"progeny") ||
       !strcmp(migstage,"both"))) {
    fprintf(stderr, "Error: \"%s\" is no valid stage of migration\n", migstage);
    exit(1);
  }
  if(!(!strcmp(multadd,"mult") || !strcmp(multadd,"add"))) {
    fprintf(stderr, "Error: \"%s\" is no valid stage of fitness accumulation\n", multadd);
    exit(1);
  }
  if(!strcmp(multadd,"add")) {
    fprintf(stderr, "Sorry, additive model isn't implemented yet. Try it with \"mult\".\n");
    exit(1);
  }

  if(staingen > finalgen){
    fprintf(stderr, "Warning: \"staining\" generation larger than final generation\n");
  }
  if(outgen > finalgen){
    fprintf(stderr, "Warning: output generation larger than final generation\n");
  }
  if (finalgen % outgen){
    fprintf(stderr, "Warning: final generation (%d) is not a multiple of\n", finalgen);
    fprintf(stderr, "         the printing steps (%d)\n", outgen);
  }
  if (!(dispersalmode==0 || dispersalmode==1)){
    fprintf(stderr,  "Error: %d is an invalid dispersal mode\n", dispersalmode);
    exit(1);
  }
  if (dispersalmode==1 && staingen > 0){
    fprintf(stderr, "Sorry, \"staining\" doesn't work for dispersalmode = 1\n");
    fprintf(stderr, "Simulation will be continued without staining.\n");
    staingen = 0;
  }
  if (!(recessivity==0 || recessivity==1 || recessivity==2 ||recessivity==3)){
    fprintf(stderr,  "Error: %d is an invalid recessivity mode for DMIs\n", recessivity);
    exit(1);
  }
  if (!(selstats==0 || selstats==1 )){
    fprintf(stderr,  "Error: %d is an invalid option for selected loci stats\n", selstats);
    exit(1);
  }
  if (selstrength < 0 || selstrength > 1){
    fprintf(stderr, "Error: selection coefficient -c %f not between 0.0 and 1.0\n", selstrength);
    exit(1);
  }
  if (ridge < 0 || ridge > 1){
    fprintf(stderr, "Error: selection along DMI ridge -a %f not between 0.0 and 1.0\n", ridge);
    exit(1);
  }

  /* ------------------- */
  /* reading input files */
  /* ------------------- */

  /* demefile */
  if (strcmp(mydemefile,"undefined") != 0){
    getdemes(mydemefile, &ndemes, &aux);
  }
  else {
    ndemes = (dispersalmode==0) ? 3 : 5; /* 3 or 5 demes as default */
    aux.cap_ad = gsl_vector_int_calloc(ndemes);
    aux.cap_pr = gsl_vector_int_calloc(ndemes);
    gsl_vector_int_set_all(aux.cap_ad, 50);
    gsl_vector_int_set_all(aux.cap_pr, 100);
  }


  /* migration matrix default */
  setmigmat(ndemes, &aux, mig);


  /* epifile */
  if (strcmp(myepifile,"undefined") != 0){
    getepi(myepifile, &nepiloci, &epi, selstrength, ridge, recessivity);
  }
  else{
    nepiloci = 0;
  }

  /* selfile */
  if (strcmp(myselfile,"undefined") != 0){
    nenviloci = getsel(myselfile,&nselloci, &sel);
  }
  else{
    nselloci = 0;
    nenviloci = 0;
  }

  nallselloci = nepiloci + nselloci;

  /* envifile */
  if (strcmp(myenvifile,"undefined") != 0  &&  nenviloci > 0 ){ /* defined and needed */
    getenvi(myenvifile, ndemes, nenviloci, &aux);
  }
  else if (strcmp(myenvifile,"undefined") == 0  &&  nenviloci > 0 ){
    /* undefined but needed, make default */
    defaultenvi(selstrength, ndemes, nenviloci, &aux);
  }
  else { /* defined or undefined but not needed */
    ;
  }


  /* ------------------- */
  /* set up output files */
  /* ------------------- */
  /* main output file */
  /* print location of markers and header to top line of file */
  fp = fopen(mainfile, "w"); /* will overwrite old file, if exist */
  if ( !fp ){
    fprintf(stderr, "Can't open %s for output!\n", mainfile);
    exit(1);
  }
  /* print header */
  fprintf(fp, "##Marker positions: "); /* line starting with #
			    will be skipped in R (read in with
			   read.table("testing.txt",sep=",",header=T)) */
  for(i=0; i<markers_per_chr; i++){
    fprintf(fp, "%.3f ", gsl_vector_get(markerloc, i));
  }
  fprintf(fp, "\nrep,gen,deme,ind,q,het,fit,njunct");
  for (j=1;j<khomologsN+1;j++){
    for (i=1;i<markers_per_chr+1;i++){
      fprintf(fp, ",l%d.%d",j,i);
    }
  }
  fprintf(fp, "\n");
  fclose(fp);

  /* haplotype output file */
  /* print location of markers and header to top line of file */
  fp = fopen(haplofile, "w"); /* will overwrite old file, if exist */
  if ( !fp ){
    fprintf(stderr, "Can't open %s for output!\n", haplofile);
    exit(1);
  }
  /* print header */
  fprintf(fp, "##Marker positions: "); /* line starting with #
			    will be skipped in R (read in with
			   read.table("testing.txt",sep=",",header=T)) */
  for(i=0; i<markers_per_chr; i++){
    fprintf(fp, "%.3f ", gsl_vector_get(markerloc, i));
  }
  fprintf(fp, "\nrep,gen,deme,ind,q,het,fit,njunct");
  for (j=1;j<khomologsN+1;j++){
    for (i=1;i<markers_per_chr+1;i++){
      fprintf(fp, ",l%d.%d",j,i);
    }
  }
  fprintf(fp, "\n");
  fclose(fp);


  /* stats output file */
  fp = fopen(statsfile, "w"); /* will overwrite old file, if exist */
  if ( !fp ){
    fprintf(stderr, "Can't open %s for output!\n", statsfile);
    exit(1);
  }
  /* print header */
  fprintf(fp, "rep,gen,deme,nind,qmean,qsd,hetmean,hetsd,fitnessmean,fitnesssd,njunctmean,njunctsd,fisdeme,fis,fst\n");
  fclose(fp);


  /* selected loci stats output file */
  if (nallselloci > 0 && selstats==1){
    fp = fopen(selstatsfile, "w"); /* will overwrite old file, if exist */
    if ( !fp ){
      fprintf(stderr, "Can't open %s for output!\n", selstatsfile);
      exit(1);
    }
    /* print header */
    fprintf(fp, "rep,gen,deme");
    for (i=0;i<nallselloci;i++){
      fprintf(fp, ",af_%d,het_%d", i,i);
    }
    fprintf(fp, "\n");
    fclose(fp);
  }
  else{
    remove(selstatsfile); /* will remove old file, if exist */
  }


  /* selected loci output file */
  if (nallselloci > 0){
    /* print header with locus info */
    fp = fopen(seloutfile, "w"); /* will overwrite old file, if exist */
    if ( !fp ){
      fprintf(stderr, "Can't open %s for output\n", seloutfile);
      exit(1);
    }
    print_selfilehead(fp, &epi, &sel, nepiloci, nselloci);
    fclose(fp);
  }
  else{
    remove(seloutfile); /* will remove old file, if exist */
  }

  /* stained allele frequencies */
  if (staingen>0) {
    /* print location of markers and header to top line of file */
    fp = fopen(stainfile, "w"); /* will overwrite old file, if exist */
    if ( !fp ){
      fprintf(stderr, "Can't open %s for output!\n", stainfile);
      exit(1);
    }
    /* print header */
    fprintf(fp, "##Marker positions: "); /* line starting with #
					    will be skipped in R (read in with
					    read.table("testing.txt",sep=",",header=T)) */
    for(i=0; i<markers_per_chr; i++){
      fprintf(fp, "%.3f ", gsl_vector_get(markerloc, i));
    }
    fprintf(fp, "\nrep,gen,deme,dye");
    for (j=1;j<khomologsN+1;j++){
      for (i=1;i<markers_per_chr+1;i++){
	fprintf(fp, ",l%d.%d",j,i);
      }
    }
    fprintf(fp, "\n");
    fclose(fp);
  }
  else{
    remove(stainfile); /* will remove old file, if exist */
  }



  /* ------------------------- */
  /* allocate memory for demes */
  /* ------------------------- */
  grid = malloc(ndemes * sizeof(struct deme));
  if(grid == NULL){
    printf("Error in allocating memory for grid (%ld, %d)\n",
	   sizeof(struct deme), ndemes);
  }



  /***************************************************************************/
  /******************************   MAIN LOOP   ******************************/
  /***************************************************************************/

  for(rep=0; rep<repN; rep++){
    /* -------------------- */
    /* initialize replicate */
    /* -------------------- */

    rng_seed = rand();
    /* rng_seed = 123; */  /* use same value to have identical stream of random numbers */
    gsl_rng_set(r, rng_seed); /* seed gsl_rng with output of rand(),
				 which was seeded with result of
				 time(NULL) */

    stain = 0;

    /* initialize each deme */
    /* (allocate memory and fill with parental individuals) */
    init_deme(grid, &aux, nenviloci);

    summarize_cur_gen(grid, nallselloci, nepiloci, nselloci, selstrength, &epi, &sel);

    /* ------------ */
    /* write output */
    /* ------------ */
    /* print generation 0 */
    gen = 0;
    fp = fopen(statsfile, "a"); /* appends */
    if ( !fp ){
      fprintf(stderr, "Can't open %s for output!\n", statsfile);
      exit(1);
    }
    write_stats(fp, rep, gen, grid);
    fclose(fp);


    /* only when printing every generation */
    if(outgen == 1){
      fp = fopen(mainfile, "a"); /* appends */
      if ( !fp ){
	fprintf(stderr, "Can't open %s for output!\n", mainfile);
	exit(1);
      }
      print_file(fp, rep, gen, markerloc, markers_per_chr, grid);
      fclose(fp);

      fp = fopen(haplofile, "a"); /* appends */
      if ( !fp ){
	fprintf(stderr, "Can't open %s for output!\n", haplofile);
	exit(1);
      }
      /* print_haplofile(fp, rep, gen, markerloc, markers_per_chr, grid); */ /* print ints */
      print_haplofile_char(fp, rep, gen, markerloc, markers_per_chr, grid); /* print chars */
      fclose(fp);

      if (nallselloci > 0){
	fp = fopen(seloutfile, "a"); /* appends */
	if ( !fp ){
	  fprintf(stderr, "Can't open %s for output!\n", seloutfile);
	  exit(1);
	}
	print_selfile(fp, grid, rep, gen,
		      &epi, &sel, nepiloci, nselloci, nallselloci);
	fclose(fp);
      }
    }

    /* selstatsfile */
    if (nallselloci > 0 && selstats==1){
      fp = fopen(selstatsfile, "a"); /* appends */
      if ( !fp ){
	fprintf(stderr, "Can't open %s for output!\n", selstatsfile);
	exit(1);
      }
      write_selstats(fp, rep, gen, grid, &epi, &sel, nepiloci, nselloci, nallselloci);
      fclose(fp);
    }


    /* ------------------------------------------ */
    /* continue for a fixed number of generations */
    /* ------------------------------------------ */

    for(gen=1; gen<=finalgen; gen++){

      if(staingen == gen){
	stain = 1; /* start with staining */
      }

      /* dummyHead_sub_next_gen.nextPtr = NULL is end of list (tail) */
      /* in new generation; chromosomes will be inserted between head (dummyHead) */
      /* and tail (NULL pointer) */
      for (k=0; k<ndemes; k++){
	grid[k].dummyHead_sub_next_gen->nextPtr = NULL;
      }


      /* make linked list (with dyn_chromosome) with all chromosomes that will be */
      /* present in next generation, after migration, recombination and selection */
      mating(grid, &aux, selstage, migstage, nallselloci,
	     &epi, &sel, nepiloci, nselloci, selstrength);

      free_cur_gen(grid); /* free memory */

      copynext2cur(grid); /* copy from linked list into chromosome vector, and clean up */

      dynamic_grid_free(grid); /* remove all elements from linked list */

      summarize_cur_gen(grid, nallselloci, nepiloci, nselloci, selstrength, &epi, &sel);


      /* ------------ */
      /* write output */
      /* ------------ */
      /* print stats every generation */
      fp = fopen(statsfile, "a"); /* appends */
      if ( !fp ){
	fprintf(stderr, "Can't open %s for output!\n", statsfile);
	exit(1);
      }
      write_stats(fp, rep, gen, grid);
      fclose(fp);


      /* selstatsfile */
      if (nallselloci > 0 && selstats==1){
	fp = fopen(selstatsfile, "a"); /* appends */
	if ( !fp ){
	  fprintf(stderr, "Can't open %s for output!\n", selstatsfile);
	  exit(1);
	}
	write_selstats(fp, rep, gen, grid, &epi, &sel, nepiloci, nselloci, nallselloci);
	fclose(fp);
      }


      /* print cur_generation to file every outgen steps */
      if (gen % outgen == 0){
	/* ind stats and genotype */
      	fp = fopen(mainfile, "a"); /* appends */
      	if ( !fp ){
      	  fprintf(stderr, "Can't open %s for output!\n", mainfile);
      	  exit(1);
      	}
     	print_file(fp, rep, gen, markerloc, markers_per_chr, grid);
      	fclose(fp);

	/* haplotype */
      	fp = fopen(haplofile, "a"); /* appends */
      	if ( !fp ){
      	  fprintf(stderr, "Can't open %s for output!\n", haplofile);
      	  exit(1);
      	}
      	/* print_haplofile(fp, rep, gen, markerloc, markers_per_chr, grid); */ /* print ints */
      	print_haplofile_char(fp, rep, gen, markerloc, markers_per_chr, grid); /* print chars */
      	fclose(fp);

	/* selected loci */
	if (nallselloci > 0){
	  fp = fopen(seloutfile, "a"); /* appends */
	  if ( !fp ){
	    fprintf(stderr, "Can't open %s for output!\n", seloutfile);
	    exit(1);
	  }
	  print_selfile(fp, grid, rep, gen,
			&epi, &sel, nepiloci, nselloci, nallselloci);
	  fclose(fp);
	}

	/* stained allele frequencies */
	if (staingen>0 && stain == 1){
	  fp = fopen(stainfile, "a"); /* appends */
	  if ( !fp ){
	    fprintf(stderr, "Can't open %s for output!\n", stainfile);
	    exit(1);
	  }
	  print_stainfile(fp, rep, gen, markerloc, markers_per_chr, grid);
	  fclose(fp);
	}

      }

    }
    free_cur_gen(grid);

  }

  /******************************************************************************/
  /****************************   END OF MAIN LOOP   ****************************/
  /******************************************************************************/

  /* --------------------------------- */
  /* write settings to info (log) file */
  /* --------------------------------- */
  char buff[20];
  struct tm *sTm;
  time_t now = time (0);
  sTm = localtime (&now);
  strftime (buff, sizeof(buff), "%Y-%b-%d %H:%M:%S", sTm);

  fp = fopen(logfile, "w"); /* will overwrite old file, if exist */
  if ( !fp ){
    fprintf(stderr, "Can't open %s for output!\n", mainfile);
    exit(1);
  }

  /* input files */
  fprintf(fp, "Run %s\n", buff);
  fprintf(fp, "Input files: deme settings '%s', epistasis '%s', single-locus '%s', environment '%s'\n", mydemefile, myepifile, myselfile, myenvifile);
  /* output files */
  fprintf(fp, "Output files: '%s'\n", mainfile);
  /* settings */
  fprintf(fp, "Number of replicates: %d\n", repN);
  fprintf(fp, "Number of generations: %d\n", finalgen);
  fprintf(fp, "Writing output every %d %s\n", outgen, ((outgen>1)?"generations":"generation"));
  fprintf(fp, "Mean fecundity: %.1f\n", (double) MeanFecund);
  fprintf(fp, "Number of homolog chromosomes: %d\n", khomologsN);
  fprintf(fp, "Migration rate: %f\n",mig);
  fprintf(fp, "Selection coefficient: %f\n", selstrength);
  fprintf(fp, "Selection along DMI ridge: %f\n", ridge);
  fprintf(fp, "Recessivity for DMIs: %d\n", recessivity);
  fprintf(fp, "Selection stage: %s\n", selstage);
  fprintf(fp, "Migration stage: %s\n", migstage);
  fprintf(fp, "Dispersal mode: %d\n", dispersalmode);
  fprintf(fp, "Markers per chromosome: %d\n", markers_per_chr);
  fprintf(fp, "Fitness model: %s\n", multadd);
  (staingen>0)?fprintf(fp, "Staining: starts at generation %d\n", staingen):
    fprintf(fp, "Staining: no\n");
  fprintf(fp, "Number of demes: %d\n", (dispersalmode==1)?ndemes-2:ndemes);
  fprintf(fp, "Adult capacity: ");
  for (k=0;k<ndemes;k++){
    fprintf(fp, "%d ",gsl_vector_int_get(aux.cap_ad,k));
  }
  fprintf(fp, "\nProgeny capacity: ");
  for (k=0;k<ndemes;k++){
    fprintf(fp, "%d ",gsl_vector_int_get(aux.cap_pr,k));
  }
  fprintf(fp, "\nMigration matrix:\n");
  for (k=0; k< ndemes; k++){ /* recipient */
    for (j=0; j< (ndemes + 2); j++){ /* source */
      fprintf(fp, "%8.3f", gsl_matrix_get(aux.migmat,k,j));
      if (j==ndemes+1){
	fprintf(fp, "\n");
      }
    }
  }
  /* marker positions */
  fprintf(fp, "Marker positions:\n");
  for(i=0; i<markers_per_chr; i++){
    fprintf(fp, "%.3f ", gsl_vector_get(markerloc, i));
  }
  fprintf(fp, "\n");

  /* some info on selected loci and environment */
  if(nallselloci == 0){
    fprintf(fp, "No selection.\n");
  }
  else{
    fprintf(fp, "Selection:\n");
    fprintf(fp, "%d epistatic loci, %d complex%s, %d-plets\n", nepiloci,
	    (nepiloci>0)?epi.ncomplexes:0, (epi.ncomplexes==1)?"":"es", (nepiloci>0)?epi.xplets:0);
    fprintf(fp, "%d selected single-loc%s, including\n", nselloci, (nselloci == 1)?"us":"i");
    fprintf(fp, "%d environmental selected loc%s\n", nenviloci, (nenviloci == 1)?"us":"i");
    if(nenviloci > 0){
      fprintf(fp, "Environment:\n");
      for(i=0; i<nenviloci; i++){
	fprintf(fp, "locus %d: ",i+1);
	for (k=0; k< ndemes; k++){
	  fprintf(fp, "%7.2f", gsl_vector_get(grid[k].environment, i));
	}
	fprintf(fp, "\n");
      }
    }
  }

  fclose(fp);

  /* ----------- */
  /* free memory */
  /* ----------- */

  gsl_vector_free(markerloc);
  gsl_rng_free(r);

  /* free grid */
  for (k=0; k<ndemes; k++){
    gsl_vector_free(grid[k].fitness);
    gsl_matrix_free(grid[k].stats);
    gsl_vector_free(grid[k].junctions);
    if(nenviloci > 0){
      gsl_vector_free(grid[k].environment);
      gsl_vector_free(grid[k].envselstrength);
    }
  }
  free(grid);

  /* other structs */
  gsl_matrix_free(aux.migmat);
  gsl_vector_int_free(aux.cap_ad);
  gsl_vector_int_free(aux.cap_pr);
  if(nenviloci > 0){
    gsl_matrix_free(aux.envi);
  }

  /* free epiloci, selloci (if they were allocated) */
  if (strcmp(myepifile,"undefined") != 0){
    gsl_vector_int_free(epi.locuscount);
    gsl_vector_int_free(epi.complex);
    gsl_vector_int_free(epi.compcount);
    gsl_vector_char_free(epi.type);
    gsl_vector_int_free(epi.chr);
    gsl_vector_free(epi.posdecimal);
    gsl_vector_char_free(epi.ancallele);
  }
  if (strcmp(myselfile,"undefined") != 0){
    gsl_vector_int_free(sel.locuscount);
    gsl_vector_char_free(sel.type);
    gsl_vector_int_free(sel.chr);
    gsl_vector_free(sel.posdecimal);
    gsl_vector_free(sel.valaa);
    gsl_vector_free(sel.valab);
    gsl_vector_free(sel.valbb);
  }


  /* finished statement */
  fprintf(stderr,"Finished run %s.\n\n", buff);

  return 0;
}





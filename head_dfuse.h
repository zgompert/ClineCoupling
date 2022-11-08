/* File: head_dfuse.h */

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
#include <gsl/gsl_rng.h>
#include <gsl/gsl_block.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>


/************************************/
/* Definition of macros/constants   */
/************************************/
/* ZG thinks khomologsN is number of chromosome, set to 1 instead of 10 */
#define khomologsN   1 
#define kstrideInd   (khomologsN * 2)
/* individuals	consist of khomologsN * 2 consecutive homologs        */
/* khomologsN * 2 is the stride to move between individuals           */

#define kchiasmaN    1  /* parameter of Poisson */
#define MeanFecund   5
/* #define MaxMig       0.5 /\* max allowed proportion of migrants for progeny dispersal */
/* 			  (not implemented yet)*\/ */


/************************************/
/* Definition of global variables   */
/************************************/

extern gsl_rng * r;  /* global random number generator */
extern int ndemes; /* number of demes */
extern int stain; /* set to true at generation x when staining starts */
extern int dispersalmode; /*  0=infinite source (default) or 1=finite source */

/************************************/
/* Definition of data structures    */
/************************************/

struct chromosome {
  int nblocks;
  gsl_vector_char * allele;
  gsl_vector_float * mapposition;
};

struct dyn_chromosome {
  int nblocks;
  gsl_vector_char * allele;
  gsl_vector_float * mapposition;
  struct dyn_chromosome * nextPtr; /* self-referential structure */
};


struct deme{
  unsigned int cap_ad; /* capacity adults */
  unsigned int cap_pr; /* capacity progeny */
  unsigned int gChromoNpop; /* all possible chromosomes in deme */
  unsigned int N; /* realized pop size */
  unsigned int Npr; /* realized number of progeny */
  gsl_vector * fitness; /* fitness for all possible inds in deme */
  gsl_matrix * stats; /* hybrid index, heterozygosity, for all possible inds in deme */
  gsl_vector * junctions; /* number of junctions per ind */
  struct chromosome * sub_cur_generation;
  struct dyn_chromosome * dummyHead_sub_next_gen;
  gsl_vector * environment; /* optimum for each environmental selected locus */
  gsl_vector * envselstrength; /* selcoeff for each environmental selected locus */
};

struct auxcont{
  gsl_vector_int * cap_ad;
  gsl_vector_int * cap_pr;
  gsl_matrix * migmat; /* migration between demes */
  gsl_matrix * envi; /* environmental selection params for demes */
};


struct epiloci{
  int ncomplexes;
  int xplets;
  gsl_vector_int *locuscount;
  gsl_vector_int *complex;
  gsl_vector_int *compcount;
  gsl_vector_char *type;
  gsl_vector_int *chr;
  gsl_vector *posdecimal;
  gsl_vector_char *ancallele;
  double dmimat[3][3]; /* defines selcoeffs for locus pairs */
  double pathmat[3][3]; /* defines selcoeffs for locus pairs */
};

struct selloci{
  gsl_vector_int *locuscount;
  gsl_vector_char *type;
  gsl_vector_int *chr;
  gsl_vector *posdecimal;
  gsl_vector *valaa;
  gsl_vector *valab;
  gsl_vector *valbb;
};



/************************************/
/* Function Prototypes 	   func.c   */
/************************************/

void usage(char * );
void getdemes(char file[], int * ndemes, struct auxcont * );
void getepi(char file[], int * loci, struct epiloci * ,
	    const double sel, const double ridge, const int recessivity);
int getsel(char file[], int * loci, struct selloci *);
void getenvi(char file[], const int ndemes, const int nenviloci, struct auxcont * );
void defaultenvi(const double selstrength,
		 const int ndemes, const int nenviloci, struct auxcont * );
void setenvi();
void setmigmat(const int ndemes, struct auxcont * , const double mig);
void setfitmat(const double sel, double dmi[][3], double path[][3], const double ridge,
	       const int recessivity);
void init_deme(struct deme * , struct auxcont * , const int nenviloci);
void pop_init(struct chromosome * cur_generation, unsigned int popbreak,
	      unsigned int gChromoNpop);
void recombine (const struct chromosome * ,
		const struct chromosome * ,
		struct chromosome * );
int init_chromosome (struct chromosome *, int);
int init_dyn_chromosome (struct dyn_chromosome * chromo, int nblocks);
int find_position (double, const struct chromosome *);
int fillgamete(gsl_vector_char *, gsl_vector_float *,
	       const struct chromosome *, const int, const int, int *);
int storehomolog (struct dyn_chromosome * homolog, struct chromosome * gamete);
int storehomologViab (struct chromosome * , struct chromosome *);
void mating (struct deme * ,
	     struct auxcont * , char * selstage, char * migstage,
	     const int nallselloci,
	     const struct epiloci * epi, const struct selloci * sel,
	     const int nepiloci, const int nselloci, const double selstrength);
void placezygote(struct chromosome gamete, struct dyn_chromosome * dummyHead);
void copynext2cur(struct deme * );
double calc_heterozygosity(const struct chromosome *,
			   const struct chromosome *);
void copy_homolog (struct chromosome *, const struct chromosome *);
int appendChiasma (gsl_vector_char * , const int , char , gsl_vector_float * , double);
void free_cur_gen(struct deme * );
void dynamic_grid_free(struct deme * );
void summarize_cur_gen(struct deme * , const int nallselloci,
		       const int nepiloci, const int nselloci, const double selstrength,
		       const struct epiloci * epi, const struct selloci * sel);
void summarize_ind(struct chromosome * ind, int, double *, double *, double *);
int getgenotype(struct chromosome * cur_generation, int, int, double);
void getselgeno(gsl_vector_int * geno, struct chromosome * cur_generation, const int ind,
		const struct epiloci * epi, const struct selloci * sel, const int nallselloci,
		const int nepiloci, const int nselloci);
int gethaplotype(struct chromosome * cur_generation, int, int, int, double);
char gethaplotype_char(struct chromosome * cur_generation, int, int, int, double);
double getfitness(gsl_vector_int * genotype, const struct epiloci * ,
		  const int nepiloci,
		  const struct selloci * , const int nselloci, const int nallselloci,
		  const double selstrength,
		  const gsl_vector * envi, const gsl_vector * envsel);
int make_parental_gamete(struct chromosome * , const char state);
void print_file(FILE * file,
		int rep, int gen, gsl_vector * markerloc, int markers_per_chr, struct deme * );
void print_haplofile(FILE * file, int rep, int gen, gsl_vector * markerloc,
		     int markers_per_chr, struct deme * );
void print_haplofile_char(FILE * file, int rep, int gen, gsl_vector * markerloc,
			  int markers_per_chr, struct deme * );
double getdyefreq(struct chromosome * cur_generation, const int N, const int chr,
		  const double position, const char dye);
void write_stats(FILE * file, int rep, int gen, struct deme * );
void write_selstats(FILE * file, const int rep, const int gen, struct deme * ,
		    const struct epiloci * epi, const struct selloci * sel,
		    const int nepiloci, const int nselloci, const int nallselloci);
void print_stainfile(FILE * file,
		     int rep, int gen, gsl_vector * markerloc, int markers_per_chr,
		     struct deme * grid);
void print_selfilehead(FILE * file,
		       const struct epiloci * epi, const struct selloci * sel,
		       const int nepiloci, const int nselloci);
void print_selfile(FILE * file, struct deme * grid, const int rep, const int gen,
		   const struct epiloci * epi, const struct selloci * sel,
		   const int nepiloci, const int nselloci, const int nallselloci);
void license(char * );

# ClineCoupling
Simulations and analyses of coupling for CSH Speciation chapter

# Some notes and theoretical background

[Barton 1983](https://www.jstor.org/stable/2408260?origin=crossref#metadata_info_tab_contents) defines the coupling coefficient as $\theta = s/r$, where $s$ is the selection coefficient and $r$ is the recombination rate between neighboring loci; $R = (L-1)/r$, where $L$ is the number of loci and $R$ is the total map length. Fitness of a hybrid heterozygous at $n$ underdominant loci is $(1-s)^n$. He notes a sharp transition between uncoupled and coupled systems (that is where $s*$, the effecitve selection experienced by a locus, approaches $S = sL$) at $\theta = 1$ when holding $S$ and $R$ constant but increasing the number of loci $L$. The divide between uncoupled and coupled is expected to be especially sharp for very large $L$. 

[Kruuk 1999](https://academic.oup.com/genetics/article/153/4/1959/6035068) revisits coupling in the context of endogenous vs exogenous selection. He focuses more on the summed coupling coefficient, defined as $\phi = (L-1)s/r$. This better reflects $s*$, that is the selection experienced by each locus, which increases with $\theta$ and with the number of loci, $L$. Kruuk notes a smooth increase in systems going from uncoupled to coupled when the per locus selection, $s$, is increased. This does not necessarily contradict the points above as the transition could be smooth as a function of $s$ and abrupt as a function of $L$. The latter, an abrupt transition as a function of $L$, is consistent with genome-wide congealing being abrupt in time as the number of selected mutations increases, e.g., [Nosil 2017](https://www.nature.com/articles/s41559-016-0001).

Our plan is to build on this past work by considering hybrid zones involving a range of values of $L$ and $\theta$ (in full factorial combinations) and specifically asking how theta and phi correspond with the variance in cline widths and centers, as well as the mean cline width, in terms of both geoegraphic and genomic clines. A primary motiviation is determing whether there is a sharp transition between uncoupled states with high cline variance and coupled states with low cline variance. Our focus on cline variance (and mean cline width) represent our interest in a parameter that can be directly quantified from genomic analyses of hybrid zones alone.

# Simulations

We are using [dfuse](https://www.uwyo.edu/buerkle/software/dfuse/) to simulation hybrid zones. The model was described in [Lindtke 2015](https://onlinelibrary.wiley.com/doi/10.1111/evo.12725). I modified the source code to allow for simple underdominance and to have only a single chromosome (see [main_dfuse.c](main_dfuse.c), [func_dfuse.c](func_dfuse.c), and [head_dfuse.h](head_dfuse.h)). This software runs individual-based simulations of secondary contact using a stepping stone model and tracks ancestry junctions.

Constant conditions for the simulations are as follows: number of demes = 110, number of generations = 2000 (log every 500), number of chromsomes = 1, migration rate between neighboring demes = 0.1, selection type = underdominance, number of genetic markers = 51.

I am varying theta and the number of loci (in all pairwise combinations). For $\theta$, I used 0.1, 0.5, 0.9, 1, 1.1, 1.5, and 2. For $L$, I used 2, 10, 100, 200, 500, and 1000. Selected loci were placed evenly along the chromosome starting at the end (i.e., there is a selected locus at each end of the chromosome). This gives $r = 1/(L-1)$. I then used $\theta = s/r$ to compute the per locus s for each simulation (the same for all loci). I ran 10 replicates of each set of simulation conditions = 7 * 6 * 10 = 300 simulations.

Here is the main submission script (in /uufs/chpc.utah.edu/common/home/gompert-group2/projects/coupling_sims/):

```bash
#!/bin/sh 
#SBATCH --time=240:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=12
#SBATCH --account=gompert-np
#SBATCH --partition=gompert-np
#SBATCH --job-name=dfuse
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=zach.gompert@usu.edu

module load gsl

cd /uufs/chpc.utah.edu/common/home/gompert-group2/projects/coupling_sims

perl fork_dfuse_main.pl unsel_*
```
Which runs:

```perl
#!/usr/bin/perl
#
#

use Parallel::ForkManager;
my $max = 80;
my $pm = Parallel::ForkManager->new($max);

@theta = (0.1,0.5,0.9,1,1.1,1.5,2); ## theta = s/r = coupling coefficient, not summed coupling coefficient

foreach $Lfile (@ARGV){ ## unsel files
	$Lfile =~ m/_(\d+)/ or die "failed to match number of Loci = L\n";
	$L = $1;
	foreach $th (@theta){
		$r = 1 / ($L-1); ## $r between neighboring loci
		$s = $th * $r; ## theta = s/r, thus s = theta * r
		foreach $j (0..9){ ## reps
			sleep 2;	
			$pm->start and next;
			$out = "o_d110_L$L"."_theta$th"."_rep$j";
			system "~/bin/dfuse_src/dfuse -d demefile_110 -s $Lfile -o $out -g 2000 -c $s -G 500 -m 0.1\n";

			$pm->finish;
		}
	}
}

$pm->wait_all_children;
```
After starting this, I decided to run a second set with m = 0.2 (that might seem high, but this is between neighboring demes and there are many demes spanning the whole space so the flux of genes across the whole system shouldn't be too crazy). I am not so much interested in different migration rates per se as this should just change the scale of everything, but I want to verify that that is true for geographic clines and verify that there really isn't much of an effect in terms of variance in clines for genomic clines.

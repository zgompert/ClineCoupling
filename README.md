# ClineCoupling
Simulations and analyses of coupling for CSH Speciation chapter

# Some notes and theoretical background

[Barton 1983](https://www.jstor.org/stable/2408260?origin=crossref#metadata_info_tab_contents) defines the coupling coefficient as $\theta = s/r$, where $s$ is the selection coefficient and $r$ is the recombination rate between neighboring loci; $R = (L-1)/r$, where $L$ is the number of loci and $R$ is the total map length. Fitness of a hybrid heterozygous at $n$ underdominant loci is $(1-s)^n$. He notes a sharp transition between uncoupled and coupled systems (that is where $s*$, the effecitve selection experienced by a locus, approaches $S = sL$) at $\theta = 1$ when holding $S$ and $R$ constant but increasing the number of loci $L$. The divide between uncoupled and coupled is expected to be especially sharp for very large $L$. 

[Kruuk 1999](https://academic.oup.com/genetics/article/153/4/1959/6035068) revisits coupling in the context of endogenous vs exogenous selection. He focuses more on the summed coupling coefficient, defined as $\phi = (L-1)s/r$. This better reflects $s*$, that is the selection experienced by each locus, which increases with $\theta$ and with the number of loci, $L$. Kruuk notes a smooth increase in systems going from uncoupled to coupled when the per locus selection, $s$, is increased. This does not necessarily contradict the points above as the transition could be smooth as a function of $s$ and abrupt as a function of $L$. The latter, an abrupt transition as a function of $L$, is consistent with genome-wide congealing being abrupt in time as the number of selected mutations increases, e.g., [Nosil 2017](https://www.nature.com/articles/s41559-016-0001).

Our plan is to build on this past work by considering hybrid zones involving a range of values of $L$ and $\theta$ (in full factorial combinations) and specifically asking how $\theta$ and $\phi$ correspond with the variance in cline widths and centers, as well as the mean cline width, in terms of both geoegraphic and genomic clines. A primary motiviation is determing whether there is a sharp transition between uncoupled states with high cline variance and coupled states with low cline variance. Our focus on cline variance (and mean cline width) represent our interest in a parameter that can be directly quantified from genomic analyses of hybrid zones alone.

# Simulations

We are using [dfuse](https://www.uwyo.edu/buerkle/software/dfuse/) to simulation hybrid zones. The model was described in [Lindtke 2015](https://onlinelibrary.wiley.com/doi/10.1111/evo.12725). I modified the source code to allow for simple underdominance and to have only a single chromosome (see [main_dfuse.c](main_dfuse.c), [func_dfuse.c](func_dfuse.c), and [head_dfuse.h](head_dfuse.h)). This software runs individual-based simulations of secondary contact using a stepping stone model and tracks ancestry junctions.

Constant conditions for the simulations are as follows: number of demes = 110, number of generations = 2000 (log every 500), number of chromsomes = 1, migration rate between neighboring demes = 0.1, selection type = underdominance, number of genetic markers = 51.

I am varying theta and the number of loci (in all pairwise combinations). For $\theta$, I used 0.1, 0.5, 0.9, 1, 1.1, 1.5, and 2. For $L$, I used 2, 10, 100, 200, 500, and 1000. Selected loci were placed evenly along the chromosome starting at the end (i.e., there is a selected locus at each end of the chromosome). This gives $r = 1/(L-1)$. I then used $\theta = s/r$ to compute the per locus s for each simulation (the same for all loci). I ran 10 replicates of each set of simulation conditions = 10 * 6 * 10 = 600 simulations (*2 for migration rate, but note that some are not well defined as fitness would have to be less than 0 to get $\theta$ given $L$, so really we have 1140 simulations).

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

@theta = (0.05,0.1,0.3,0.5,0.7,0.9,1,1.1,1.5,2); ## theta = s/r = coupling coefficient, not summed coupling coefficient

foreach $Lfile (@ARGV){ ## unsel files
	$Lfile =~ m/_(\d+)/ or die "failed to match number of Loci = L\n";
	$L = $1;
	foreach $th (@theta){
		$r = 1 / ($L-1); ## $r between neighboring loci
		$s = $th * $r; ## theta = s/r, thus s = theta * r
		foreach $j (0..9){ ## reps
                        sleep 2;        
                        $pm->start and next;
                        $out = "o_m2_d110_L$L"."_theta$th"."_rep$j";
                        system "~/bin/dfuse_src/dfuse -d demefile_110 -s $Lfile -o $out -g 2000 -c $s -G 500 -m 0.2\n";
                        $out = "o_d110_L$L"."_theta$th"."_rep$j";
                        system "~/bin/dfuse_src/dfuse -d demefile_110 -s $Lfile -o $out -g 2000 -c $s -G 500 -m 0.1\n";


			$pm->finish;
		}
	}
}

$pm->wait_all_children;
```
After starting this, I decided to run a second set with m = 0.2 (incorporated in code above; that might seem high, but this is between neighboring demes and there are many demes spanning the whole space so the flux of genes across the whole system shouldn't be too crazy). I am not so much interested in different migration rates per se as this should just change the scale of everything, but I want to verify that that is true for geographic clines and verify that there really isn't much of an effect in terms of variance in clines for genomic clines.

# Simulation results: cline in hybrid indexes

As a first summary of the simulation results, I plotted geographic clines in the mean hybrid index for each site ([SummarizeHi.R](SummarizeHi.R)). The results are shown here, [clinesHi.pdf](https://github.com/zgompert/ClineCoupling/files/9984685/clinesHi.pdf). Each panel is a deparate set of simulation conditions; the title gives the number of demes (110) and $\theta$ (migration 0.1 or, with m2 = 0.2). Each line is a replicate (not differnt loci) and black versus gray denote 2000 versus 1500 generations. Key take homes are that 2000 generations appears to be sufficient to reach equilibirum (not notably different from 1500), and thus will work for downstream analyses, and that the overall shape of the hybrid index cline is consistent across replicates and looks to be closely tied to $\theta$.

# Geographic clines for simulations

The simulations have a strong geographic axis and a common scale (at least one for each migration rate) and thus first quantifying variation in geographic clines makes great sense here. That is what we did. One issue is the different mathematical forms expected for single locus (sigmoidal) vs multilocus (stepped) clines. Fortunately, in both cases the center part of the cline should reflect the total selection on the locus and should have the same form. Specifically, we can estimate the maximum gradient of the cline by tranforming $p$ (allele frequency) to a logit scale at which point the relationship between location and $\mathrm{logit}(p)$ should be approximately linear (e.g., [Kruuk 1999](https://academic.oup.com/genetics/article/153/4/1959/6035068),  [Barton 1989](https://www.nature.com/articles/341497a0)). The cline width for $p$ is related to the slope on the logit scale as $w = \frac{1}{0.25 \beta}$ Kruuk recommends estimating the gradient between populations with $p \approx 0.12$ and $p \approx 0.88$ (-2 and 2 on the logit scale). I did something a bit different with the same aim. Specifically, in some cases the transition occurs between a single pair of populations or even between a pair of populations, which renders this sometimes infeasible. Instead, I find the point where $p$ for each locus is closest to 0.5. The, I take that popualtion and five on either side as the set to estimate the maximum gradient. Thus, I always have 11 data points to fit the linear regression. This means sometimes I might include populations outside of the linear component (where the gradient is maxium) but this seems like a reasonable overall tradeoff in terms of automating the analysis of 1140 data sets.

The main R script for running these analyses is [fitGeoClines.R](fitGeoClines.R). This uses the hierarchical Bayesian model defined in [geoclinemod.stan][geoclinemod.stan](geoclinemod.stan), which we fit with Hamiltonian Monte Carlo using `rstan` (version 2.21.7, GitRev: 2e1f913d3ca3) (the script also includes a non-Bayesian version of the analysis that I played with first). This was done with `R` version 4.1.3. For each data set, I fit the cline models with 20 chains using a 1000 iteration warmup and 1200 iterations total. I then extracted the mean, standard deviation and coefficient of variation of the cline widths, and the actual hierarchical mean and standard deviation on the logit scale (gradient for logit $p$). 

As expected, the mean and variance in cline widths and gradients are strongly associated with the $\Theta$, the coupling coefficient (see the linear model fits in [fitGeoClines.R](fitGeoClines.R)). The $r^2$ for explaining $\Theta$ are on around 0.7 to 0.8. This is even better than we do for genomic clines, though I suspect this reflects being able to take advantage of the high $\Theta$ simulations here in a way that we really can't (no hybrids) for genomic clines (see the next section for details). Still, the key advantage of the genomic clines approach is that the scale (hybrid index) is the same for all hybrid zones, whereas the geographic scale depends on the scale of dispersal. In other words, a high or low SD in terms of geographic clines depends on an (often) unknown geographic scale of dispersal, which isn't true for the genomic clines. This isn't an issue for the simulations as all involve 110 demes (though the migration rate does vary from 0.1 to 0.2 between demes).

# Genomic clines for simulations

I am using the logit-logistric function for genomic clines following [Fitzpatrick 2013](https://onlinelibrary.wiley.com/doi/full/10.1002/ece3.609) (but modeling ancestry unlike the original). This avoids having to splice functions in `rstan`. The version for the simulations assumes hybrid index and parental allele frequencies are known. I have versions with and without sum-to-zero constraints on the cline parameters. I am trying without first. The model has two cline parameter $u$ and $v$ and takes the form $\Phi = \frac{h^v}{h^v + (1-h)^v * e^u}$, where $\Phi$ is the ancestry probability, and $h$ is hybrid index (indexed by individual). $v$ measures the slope (gradient) relative to genome-wide average (1 denotes same as genome-wide avaerage) and $u$ is related to the center (both are indexed by locus). To scale the center such that it is the hybrid index at which $\Phi = .5$, I am using $\mathrm{logit}(c) = \frac{u}{v}$ [Bailey 2022](https://www.authorea.com/doi/full/10.22541/au.164848698.82546348). I am placing hierarchical nomrla priors on $\mathrm{logit}(c)$ and $\mathrm{log10}(v)$ so that both have set, expected means of 0 (which is the expectation for genome genome-average admixture). Thoe main goasl is then to estimate the hierarchical prior SDs. The stan models are [simple_clinemod_sz.stan](simple_clinemod_sz.stan) and [simple_clinemod_nosz.stan](simple_clinemod_nosz.stan). My R script for fitting the cline models is [fitGenomicClines.R](fitGenomicClines.R). Here, I am keeping demes with mean hybrid index between 0.1 and 0.9 and then sampling up to 300 hybrids from those demes (I am not even attempting to fit models where fewer than 50 individuals are available, and even in these cases most might not be hybrids, but rather just in a deme with a mean hybrid index within my specified bounds). Each fit involves four Hamiltonian MC chains with 2000 iterations and a 1000 iteration warmup.
 
I am parallelizing the cline fitting with each job running a set of 10 data sets:

```perl
#!/usr/bin/perl
#
# fit genomic cline models 
#

use Parallel::ForkManager;
my $max = 10;
my $pm = Parallel::ForkManager->new($max);

foreach $j (1..114){
	$pm->start and next;
	system "R CMD BATCH --no-save --no-restore \"--args $j\" fitGenomicClines.R\n";
	$pm->finish;
}
$pm->wait_all_children;
```
Next, I combined all of the output in a single `R` workspace, visualized all of the clines, and quantified the relationship between the variability of cline center and width and the coupling coefficient $\Theta$, see [combinePlotClines.R](combinePlotClines.R). A few things to keep in mind. Clines were only fit when 50 or more individuals could be sampled from demes with a mean hybrid index between 0.1 and 0.9. Even in those cases where this was true, there were sometimes very, very few hybrids with intermediate ancestry (i.e., individuals with hybrid index between 0.1 and 0.9). Specifically, across all 1140 data sets, the mean number of hybrids by this definition was 765 (meidan 405), but 44 data sets had no hybrid by this definition and 200 data sets had fewer than 10. These all correspond to cases with large coupling coefficients ($>$1). Clines fit from these data sets (when they were fit) appear less reliable, especially in terms of the soft centering (though some with $\sim$10 hybrids look okay). This drives home the point that you need a reasonable hybrid index gradient (not a giant gap with no hybrids) to fit genomic cline modeles. Thus, for my summaries below (and see [combinePlotClines.R](combinePlotClines.R)), I only considered simulations with more than 10 or more than 50 hybrids (available total, not necessarily sampled for the analysis).

Coupling explains a notable proportion of the variation in the SD for cline center and width (as parameters of the hierarchical distribtions defined in the `stan` model), with variability in cline width better explaing coupling than cline center (see [CouplingVsClineSD.pdf](https://github.com/zgompert/ClineCoupling/files/10027328/CouplingVsClineSD.pdf)). These models include linear and quadratic effects of $\Theta$. I don't see an obvious, sharp divide in SD by $\Theta$, but hybrids clearly become rarer beyond $\Theta = 1$ or so and by $\Theta = 2$ (mostly not shown) we rarely have a sufficient number of hybrids to fit the model. The cases where we can fit the model suggest a very small SD. Flipping things, I asked how well both SDs predict the coupling coefficient as this would be the direction taken empirically. The results are reasonably promising (though far from perfect), with an $r^2$ of about 0.5. The SDs explain less of the variation in ths summed coupling coefficient, $\phi = L\Theta$, with $r^2$ of about 0.3. 

Next steps: I am mostly happy with this but could imagine a finer grid of coupling coefficients but stopping at 1.1 or so. Also, I want to summarize the prevalence of hybrids (distribution of hybrid index or something like that) as another metric for getting at coupling (i.e., when $\Theta$ is high, it is the absence of hybrids rather than clines per se that is most striking). 

# Tests of robustness

I ran and analyzed simulated data sets to assess the effect of sampling on the inferred cline SDs and thus the associated coupling coefficient estimates. My default conditions (what we used for the core simulations) involved sampling up to 300 individuals from the core hybrid zone. Here, I defined this as demes where there average hybrid index was between 0.1 and 0.9. Those are our base conditions. For the comparison, I used the same simulations but only sampled from demes on the edges of the hybrid zone, specifically those demes where the mean hybrid index was between 0.1 and 0.2 or 0.8 and 0.9 (i.e., individuals mostly similar to one or the other parent). I then compared the cline SD estimates for the two sets of samples for each data set. See [fitCompareClines.R](fitCompareClines.R) for the main fits and [combineComps.R](combineComps.R) for the summary analysis. Note that I did not analyze all of the original data sets, but focused on intermediate coupling coefficients (good numbers of hybrids) so that subsampling would be possible.

The SDs for center and slope/width are highly correlated for the two samples (Pearson correlation ~0.9), but there is a general tendency for the poorer sampling (edge sampling) to result in increased SDs for center and slope (which corresponds to lower coupling coefficients). This makes sense from theory in the sense that coupling has the biggest effect of introgression in the center of a hybrid zone (the stepped part of a geographic cline) and less of an effect as you move away from the center (as LD among foreign alleles break down... where you get exponential decay in classic mulitlocus clines). Thus, by focusing on the edges one might expect a greater SD in clines as coupling is less influential. The results are summarized in [F_SamplingEffect.pdf](https://github.com/zgompert/ClineCoupling/files/10356086/F_SamplingEffect.pdf). 

Predictive power for single deme simulations...

# Summary of results

The data sets we are working with are summarized [here](https://docs.google.com/document/d/13ojSDHTLW1YxPb16ddE27iLHJs0N7P_foxr9rGC1Vnk/edit) and are on [dropbox](https://www.dropbox.com/scl/fo/axlbcw8yhpnhy8oktpgb3/h?dl=0&rlkey=eqgxopr3ynx6ylm5vv67wr7pv). Below I am keeping track of the SD parameter estimates for the data sets that we have finished analyzing. Quick note, keep an eye on how well the soft centering is going and whether it might be necessary to account for how well it is working in the SDs (looking pretty damn good so far).

| Organism | Subset of loci | $\sigma_c$ | $\sigma_v$ | SD $c$ | SD $v$ | $\theta$ | Script |
|----------|----------------|------------|------------|--------|--------|---------|--------|
| *Agalychnis* | Anc. info. and missing | 1.42 | 0.56 | 1.42 | 0.39 | 0.12 | [FitClineModel_Agalychnis.R](FitClineModel_Agalychnis.R) |
| *Aloutatta* | 1000 anc. informative | 1.02 | 0.23 | 1.01 | 0.20 | 0.65 | [FitClineModel_Alouatta.R](FitClineModel_Alouatta.R) |
| *Ceononympha* |  All anc. informative | 1.11 | 0.31 | 1.11 | 0.29 | 0.44 | [FitClineModel_Ceononympha.R](FitClineModel_Ceononympha.R) |
| *Corvus* | Anc. info. and missing | 0.86 | 0.40 | 0.86 | 0.34 | 0.14 | [FitClineModel_Corvus.R](FitClineModel_Corvus.R) |
| *Croatalus* | 1000 anc. and missing | 0.91 | 0.33 | 0.90 | 0.29 | 0.37 | [FitClineModel_Croatalus.R](FitClineModel_Croatalus.R) |
| *Encelia* | 1000 anc. and missing | 0.38 | 0.31 | 0.27 | 0.38 | 0.41 | [FitClineModel_Encelia_ASPVEN.R](FitClineModel_Encelia_ASPVEN.R) |
| *Fundulus* | 1000 anc. and missing | 0.45 | 0.12 | 0.45 | 0.11 | 1.30 | [FitClineModel_Fundulus.R](FitClineModel_Fundulus.R) | 
| *Fundulus* Coss| 1000 anc. and missing | 0.33 | 0.08 | 0.34 | 0.08 | 1.58 | [FitClineModel_Fundulus_Cossatot.R](FitClineModel_Fundulus_Cossatot.R) | 
| *Gryllus* CT | All anc. informative | 0.81 | 0.24 | 0.81 | 0.22 | 0.67 | [FitClineModel_Gryllus_CT.R](FitClineModel_Gryllus_CT.R) |
| *Gryllus* PA | All anc. informative | 0.48 | 0.15 | 0.48 | 0.14 | 1.15 | [FitClineModel_Gryllus_PA.R](FitClineModel_Gryllus_PA.R) |
| *Hirundo* H1 | All anc. informative | 0.12 | 0.39 | 0.12 | 0.36 | -0.07 | [FitClineModel_Hirundo_H1.R](FitClineModel_Hirundo_H1.R) |
| *Hirundo* H2 | All anc. informative | 0.53 | 0.31 | 0.52 | 0.30 | 0.42 | [FitClineModel_Hirundo_H2.R](FitClineModel_Hirundo_H2.R) |
| *Lissotriton* L | Anc. info. and missing | 0.35 | 0.14 | 0.35 | 0.14 | 1.26 | [FitClineModel_lissotriton_L.R](FitClineModel_lissotriton_L.R) |
| *Lissotriton* R | Anc. info. and missing | 0.97 | 0.14 | 0.97 | 0.14 | 0.92 | [FitClineModel_lissotriton_R.R](FitClineModel_lissotriton_R.R) |
| *Lycaeides* | All anc. informative | 0.94 | 0.35 | 0.93 | 0.32 | 0.31 | [FitClineModel_Lycaeides.R](FitClineModel_Lycaeides.R) |
| *Lycaeides* | Autosomal only | 0.71 | 0.32 | 0.71 | 0.29 | 0.39 | [FitClineModel_Lycaeides.R](FitClineModel_Lycaeides.R) |
| *Lycaeides* | Z only | 1.50 | 0.36 | 1.46 | 0.32 | 0.41 | [FitClineModel_Lycaeides.R](FitClineModel_Lycaeides.R) |
| *Motacilla* | 1000 anc. and missing | 0.78 | 0.49 | 0.75 | 0.40 | -0.21 | [FitClineModel_Motacilla.R](FitClineModel_Motacilla.R) |
| *Mus* BV | 1000 anc. and missing | 0.69 | 0.16 | 0.69 | 0.15 | 1.00 | [FitClineModel_Mus_BV.R](FitClineModel_Mus_BV.R) |
| *Mus* CZ | 1000 anc. and missing | 0.70 | 0.14 | 0.70 | 0.13 | 1.07 | [FitClineModel_Mus_CZ.R](FitClineModel_Mus_CZ.R) |
| *Mus* SX | 1000 anc. and missing | 1.01 | 0.23 | 1.01 | 0.19 | 0.66 | [FitClineModel_Mus_SX.R](FitClineModel_Mus_SX.R) |
| *Mytilus* | Anc. info. and missing | 0.29 | 0.27 | 0.29 | 0.22 | 0.62 | [FitClineModel_Mytilus.R](FitClineModel_Mytilus.R) |
| *Nematocharax* | Anc. info. and missing | 0.80 | 0.41 | 0.80 | 0.34 | 0.08 | [FitClineModel_Nematocharax.R](FitClineModel_Nematocharax.R) |
| *Oleria* | 1000 anc., het. and miss. | 1.73 | 0.50 | 1.72 | 0.41 | 0.42 | [FitClineModel_Oleria.R](FitClineModel_Oleria.R) |
| *Papilio* | All anc. informative | 0.29 | 0.21 | 0.29 | 0.18 | 0.93 | [FitClineModel_Papilio.R](FitClineModel_Papilio.R) | 
| *Papio* | Anc. info., miss. and het. | 0.84 | 0.41 | 0.84 | 0.36 | 0.10 | [FitClineModel_Papio.R](FitClineModel_Papio.R) | 
| *Picea* | All anc. informative | 0.67 | 0.22 | 0.67 | 0.21 | 0.77 | [FitClineModel_Picea_glauXstich.R](FitClineModel_Picea_glauXstich.R) |
| *Pinus* | Anc. info. and missing | 1.10 | 0.37 | 1.10 | 0.32 | 0.29 | [FitClineModel_Pinus.R](FitClineModel_Pinus.R) |
| *Poecile* | 1000 anc. and missing | 0.94 | 0.40 | 0.94 | 0.37 | 0.17 | [FitClineModel_Poecile.R](FitClineModel_Poecile.R) |
| *Poecile* MO | 1000 anc. and missing | 0.82 | 0.38 | 0.81 | 0.35 | 0.19 | [FitClineModel_Poecile_MO.R](FitClineModel_Poecile_MO.R) |
| *Sceloporus* | Anc. info. and missing | 1.09 | 0.21 | 1.09 | 0.21 | 0.69 | [FitClineModel_Sceloporus.R](FitClineModel_Sceloporus.R) |
| *Sternotherus* | Anc. info. and missing | 0.57 | 0.28 | 0.56 | 0.27 | 0.55 | [FitClineModel_Sternotherus.R](FitClineModel_Sternotherus.R) |
| *Zonotrichia* | Anc. info. and missing | 0.63 | 0.57 | 0.63 | 0.42 | -0.63 | [FitClineModel_Zonotrichia.R](FitClineModel_Zonotrichia.R) | 


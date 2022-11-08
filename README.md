# ClineCoupling
Simulations and analyses of coupling for CSH Speciation chapter

# Some notes and theoretical background

[Barton 1983](https://www.jstor.org/stable/2408260?origin=crossref#metadata_info_tab_contents) defines the coupling coefficient as thete = s/r, where is the selection coefficient and r is the recombination rate between neighboring loci, R = (L-1)/r, where L is the number of loci and R is the total map length. He notes a sharp transition between uncoupled and coupled systems (that is where s* approaches S) at theta = 1 when holding S and R constant but increasing the number of loci L. This is divide is expected to be especially shart for very large L. 

[Kruuk 1999](https://academic.oup.com/genetics/article/153/4/1959/6035068) revisits coupling in the context of endogenous vs exogenous selection. He focuses more on the summed coupling coefficient, defined as Phi = (L-1)s/r. This better reflects s*, that is the selection experienced by each locus, which increases with theta and with the number of loci, L. Kruuk notes a smooth increase in systems going from uncoupled to coupled when the per locus selection, s, is increased. This does not necessarily contradict the points above as the transition could be smooth as a function of s and abrupt as a function of L. The latter, an abrupt transition as a function of L, is consistent with genome-wide congealing being abrupt in time as the number of selected mutations increases, e.g., [Nosil 2017](https://www.nature.com/articles/s41559-016-0001).

Our plan is to build on this past work by considering hybrid zones involving a range of values of L and theta (in full factorial combinations) and specifically asking how theta and phi correspond with the variance in cline widths and centers, as well as the mean cline width, in terms of both geoegraphic and genomic clines. A primary motiviation is determing whether there is a sharp transition between uncoupled states with high cline variance and coupled states with low cline variance. Our focus on cline variance (and mean cline width) represent our interest in a parameter that can be directly quantified from genomic analyses of hybrid zones alone.

# Simulations

We are using [dfuse](https://www.uwyo.edu/buerkle/software/dfuse/) to simulation hybrid zones. The model was described in [Lindtke 2015](https://onlinelibrary.wiley.com/doi/10.1111/evo.12725). I modified the source code to allow for simple underdominance and to have only a single chromosome (see ). 

functions {
	real calc_phi(real h, real vv, real uu){
		real phi;
		phi = (h^vv)/((h^vv)+((1-h)^vv)*exp(uu));
		return phi;
	}
	real calc_lik(real g, real p0, real p1, real h, real vv, real uu){
		real prob;
		real phi;
		phi = calc_phi(h, vv, uu);
                if(g==0)
                        prob = log(phi * (1-p1) + (1-phi) * (1-p0)) + log(phi * (1-p1) + (1-phi) * (1-p0));
                else if (g==1) 
                        prob = log(phi * (1-p1) + (1-phi) * (1-p0)) + log(phi * p1 + (1-phi) * p0);
                else    
                        prob = log(phi * p1 + (1-phi) * p0) + log(phi * p1 + (1-phi) * p0);


		return prob;
	}
}

data{
	int L; /* # of loci */
	int N; /* # of organisms */
	real<lower=0, upper=2> G[N, L]; /* matrix of G*/
	vector<lower=0, upper=1>[N] H; /* vector of hybrid indexes */
	vector<lower=0, upper=1>[L] P0; /* parent 0 allele frequencies */
	vector<lower=0, upper=1>[L] P1; /* parent 0 allele frequencies */
}

parameters{
	vector<lower=0.001,upper=0.999>[L-1] center_raw; /* cline center parameter */
	vector<lower=0.1,upper=10>[L-1] v_raw; /* cline width parameter */
	real<lower=0> sc; /* sigma for center*/
	real<lower=0> sv; /* sigma for v*/
}

transformed parameters{
	vector[L] center = append_row(center_raw, 1/(1+exp(sum(logit(center_raw))))); /* cline center parameter */
	vector[L] v = append_row(v_raw, pow(10,-sum(log10(v_raw)))); /* cline width parameter */
	vector<lower=-100,upper=100>[L] u; /* cline u parameter */
	for(i in 1:L){
		u[i] = logit(center[i]) * v[i];
	}
}

model{

	for(i in 1:L){
		for(j in 1:N){
			/* increment likelihood */
			target += calc_lik(G[j,i], P0[i], P1[i], H[j], v[i], u[i]);
		}
	}
	for(i in 1:(L)){
		/* increment prior on v and u */
		target += normal_lpdf(logit(center[i]) | 0, sc); /* make sure 0 works for both */
		target += normal_lpdf(log10(v[i]) | 0, sv);

	}
	target += normal_lpdf(sc | 0, 1);
	target += normal_lpdf(sv | 0, 1);
}


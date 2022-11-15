data{
	int L; /* # of loci */
	int N; /* # of demes */
	real P[N, L]; /* matrix of logit allele frequencies */
}

parameters{
	vector<lower=-100,upper=100>[L] beta; /* vector of slopes*/
	vector<lower=-10,upper=10>[L] beta0; /* vector of intercepts */
	real mu; /* mean for w*/
	real<lower=0> sigma; /* sigma for w*/
	real<lower=0> serr; /* residual sd */
}

transformed parameters{

	vector[L] w; /* vector of widths */
	real wbar;
	real wsd;
	real wcv; /* coefficient of variation for w */
	for(i in 1:L){
		w[i] = 1/(0.25 * beta[i]);
	}
	wbar = mean(w);
	wsd = sd(w);
	wcv = wsd/fabs(wbar); /* coefficient of variation based on mean and sd*/
}


model{

	for(i in 1:L){
		for(j in 1:N){
			/* increment likelihood */
			target += normal_lpdf(P[j,i] | beta0[i] + beta[i] * j, serr); 
		}
	}
	for(i in 1:L){
		/* increment prior on v and u */
		target += normal_lpdf(beta[i] | mu, sigma);
		target += normal_lpdf(beta0[i] | 0, 5);
	}
	target += normal_lpdf(mu | 0, 5);
	target += gamma_lpdf(sigma | 0.1, 0.01);
	target += gamma_lpdf(serr | 0.1, 0.01);
}


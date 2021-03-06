data {
	int<lower=1> Y;          // The total number of years in the time series that can be used 
	real<lower=0> N9[Y];   // # of neophytes 9 years ago. although integers, modeled as a real
	real<lower=0> N10[Y];  // # of neophytes 10 years ago. although integers, modeled as a real
	real<lower=0> N11[Y];  // # of neophytes 11 years ago. although integers, modeled as a real
	real<lower=0> N12[Y];  // # of neophytes 12 years ago. although integers, modeled as a real
	real<lower=0> H[Y];    // # of observed hatchlings. although integers, modeled as a real
	
}

parameters {
	real<lower=0, upper = 10> sigma;    // standard deviation of log-normal hatchling counts
	//real<lower=0, upper = 1.0> phi[Y];     // annual average survival rates for juveniles
	real<lower=0, upper = 1.0> phi;     // annual average survival rates for juveniles
	//real<lower=0, upper = 1.0> p[4];    // proportions of neophytes that hatched x yeas ago, x = 9-12 years
	real<lower=0, upper=1.0> p[Y, 4];
}

transformed parameters {
	real alpha = 1.0;        // Dirichlet prior parameters
	real<lower=0> mu[Y];     // the log-mean of the # hatchlings 
	
	for (i in 1:Y){
		mu[i] = log(p[i, 1] * N9[i] * pow(phi,-9) + 
		             p[i, 2] * N10[i] * pow(phi,-10) +
		             p[i, 3] * N11[i] * pow(phi,-11) + 
		             p[i, 4] * N12[i] * pow(phi,-12)); 

		/*mu[i] = log( p[1] * N9[i] * pow(phi[i],-9) + 
		             p[2] * N10[i] * pow(phi[i],-10) +
		            p[3] * N11[i] * pow(phi[i],-11) + 
		            p[4] * N12[i] * pow(phi[i],-12)); */

		/*mu[i] = log(p[1] * N9[i] * pow(phi,-9) + 
		            p[2] * N10[i] * pow(phi,-10) +
		            p[3] * N11[i] * pow(phi,-11) + 
		            p[4] * N12[i] * pow(phi,-12));  */
	}	

	
}

model
{
	// Priors     	
	sigma ~ uniform(0 , 10);
	phi ~ beta(1,1);

	//for (i in 1:4){
	//	p[i] ~ beta(1, 1);
	//}

	for (i in 1:Y){
		//phi[i] ~ beta(1,1);
		H[i] ~ lognormal(mu[i], sigma);
		for (j in 1:4){
			p[i, j] ~ beta(1, 1);
		}

	}
  
}

generated quantities {
  // log_lik is for use with the loo package
  vector[Y] log_lik;
  for(i in 1:Y) {
  	log_lik[i] = lognormal_lpdf(H[i] | mu[i], sigma);
  }
}


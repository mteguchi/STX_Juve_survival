data {
	int<lower=1> Y;            // The total number of years in the time series that can be used 
	int<lower = 1> n_ages;             // the number of age groups in AFR
	real<lower=0> N[Y+n_ages];      // # of neophytes. although integers, modeled as a real
	real<lower=0> H[Y];        // # of observed hatchlings. although integers, modeled as a real
	
}

parameters {
	real<lower=0, upper = 10> sigma;       // standard deviation of log-normal hatchling counts
	real<lower=0, upper = 1.0> phi;        // annual average survival rates for juveniles	
	//real<lower=0, upper=1.0> p[Y, 4];
	simplex[n_ages] q;
}

transformed parameters {
	real<lower=0> mu[Y];     // the log-mean of the # hatchlings 
	real<lower=0> alpha = 1.0;   // dirichlet parameters
	//real tmp;
	
	for (i in 1:Y){
		//tmp = 0.0;

		//for (j in 1:n_ages){
		//	tmp += q[j] * N[i + j -1] * pow(phi, -(n_ages + j));		
		//}
		//mu[i] = log(tmp);

		mu[i] = log(q[1] * N[i] * pow(phi,-12) +
						q[2] * N[i+1] * pow(phi,-11) +
						q[3] * N[i+2] * pow(phi,-10) + 
		             q[4] * N[i+3] * pow(phi,-9) +
		             q[5] * N[i+4] * pow(phi,-8) + 
		             q[6] * N[i+5] * pow(phi,-7)); 

	}	

	
}

model
{
	// Priors     	
	sigma ~ uniform(0 , 10);
	phi ~ beta(1,1);
	q ~ dirichlet(rep_vector(alpha, n_ages));

	for (i in 1:Y){
		//phi[i] ~ beta(1,1);
		H[i] ~ lognormal(mu[i], sigma);		

	}
  
}

generated quantities {
  // log_lik is for use with the loo package
  vector[Y] log_lik;
  for(i in 1:Y) {
  	log_lik[i] = lognormal_lpdf(H[i] | mu[i], sigma);
  }
}


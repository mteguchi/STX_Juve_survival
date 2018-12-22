// A Stan model for estimating time-specific survival rate (phi) and time-specific
// proportions of neophytes (q) that mature at different ages. This model
// is different from Model_qt_phit.stan by having the flexible number of ages that
// are considered for ages at first reproduction (n_ages). The maximum age can be set 
// through max_age variable. 

// Neophytes vector is in reverse chronological order to accomodate the model. So...
// Estimated hatchlings (mu) are also in the reverse order. 

data {
	int<lower=1> Y;               // The total number of years in the time series that can be used 
	int<lower = 1> n_ages;             // the number of age groups in AFR
	real<lower=0> N[Y+n_ages-1];      // # of neophytes. although integers, modeled as a real
	real<lower=0> H[Y];           // # of observed hatchlings. although integers, modeled as a real
	real<lower=0> max_age;        // maximum age at first reproduction
}

parameters {
	real<lower=0, upper = 10> sigma;       // standard deviation of log-normal hatchling counts
	real<lower=0, upper = 1.0> phi[Y];        // annual average survival rates for juveniles	
	simplex[n_ages] q[Y+n_ages-1];
}

transformed parameters {
	real<lower=0> mu[Y];     // the log-mean of the # hatchlings 
	real<lower=0> alpha = 1.0;   // dirichlet parameters
	real tmp;

	for (i in 1:Y){
		tmp = 0.0;

		for (j in 1:n_ages){
			tmp += q[i + j -1, j] * N[i + j -1] * pow(phi[i], -(max_age - j + 1));		
		}
		mu[i] = log(tmp);

		/*mu[i] = log(q[i, 1] * N[i] * pow(phi[i],-9) + 
		             q[i+1, 2] * N[i+1] * pow(phi[i],-10) +
		             q[i+2, 3] * N[i+2] * pow(phi[i],-11) + 
		             q[i+3, 4] * N[i+3] * pow(phi[i],-12)); */

	}	

	
}

model
{
	// Priors     	
	sigma ~ uniform(0 , 10);
	//phi ~ beta(1,1);

	for (i in 1:Y){
		phi[i] ~ beta(1,1);
		H[i] ~ lognormal(mu[i], sigma);
		q[i] ~ dirichlet(rep_vector(alpha, n_ages));
		

	}
  
}

generated quantities {
  // log_lik is for use with the loo package
  vector[Y] log_lik;
  for(i in 1:Y) {
  	log_lik[i] = lognormal_lpdf(H[i] | mu[i], sigma);
  }
}


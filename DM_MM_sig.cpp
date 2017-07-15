#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_INTEGER(disp_model); // 1 = half-normal, 2 = exponential, 3 = half-cauchy
  DATA_INTEGER(count_model); // 1 = poisson, 2 = negative binomial
  DATA_INTEGER(is_random);
  DATA_INTEGER(is_asymptotic);
  // DATA_INTEGER(ncreeks);  // ncreeks is the number of release locations (creeks). Not currently used.
  DATA_INTEGER(nrel);  // nrel is the number of tagged fish released
  DATA_INTEGER(nsites); // the number of trapping sites at each location
  DATA_INTEGER(nperiods); // the number of recap events at each location
  DATA_INTEGER(ntraps); // the number of traps at each site
  DATA_MATRIX(countmat); //(1,nsites*nperiods,1,ntraps) // the number of recap events at each location
  DATA_VECTOR(distances); //(1,nsites) // vector of distances from release site to center of recap site
  DATA_VECTOR(times); //(1,nperiods) // vector of times since release
  DATA_SCALAR(dist_cutoff);

  PARAMETER(logit_survival);
  PARAMETER(log_detectability);
  PARAMETER(log_sig_disp_mu)
  PARAMETER(log_sig_disp_sig);
  PARAMETER_VECTOR(log_sig_disp_eps);
  PARAMETER(overdispersion);
  
  Type survival = exp(logit_survival)/(1+exp(logit_survival));
  Type detectability = exp(log_detectability);
  Type sig_disp_sig = exp(log_sig_disp_sig);
  vector<Type> sig_disp(nperiods);
  // Type overdispersion = exp(log_overdispersion);

  vector<Type> predcount(nsites * nperiods);// * ntraps);
  vector<Type> obscount(nsites * nperiods);// * ntraps);
  vector<Type> dist_factor(nsites * nperiods);// * ntraps);
  matrix<Type> outmat(nsites * nperiods, 5);// * ntraps, 5);
  vector<Type> fifty_pct(nperiods);
  vector<Type> pct_at_dist(nperiods);
    
  Type f = 0; // objective function value
  
  int counter = 0;
  for (int i=0; i<nperiods; i++) {
    f -= dnorm(log_sig_disp_eps(i), Type(0.0), sig_disp_sig, true) * is_random;
    if(is_asymptotic == 0) {
      sig_disp(i) = exp(log_sig_disp_mu + log_sig_disp_eps(i));
    } else {
      sig_disp(i) = exp(log_sig_disp_mu)*times(i)/(exp(log_sig_disp_sig)+times(i));
    }
    for (int j=0; j<nsites; j++) {
      // pool traps
      obscount(counter) = 0;
      for (int k=0; k<ntraps; k++) {
        obscount(counter) += countmat(j+(i*nsites), k);
      }
      
      // set up dispersal
      if(disp_model == 1) { // 1 = half-normal
        // dist_factor(counter) = Type(2.0) * (pnorm(distances(j)+site_width, Type(0.0), sig_disp) - 
        //     pnorm(distances(j)-site_width, Type(0.0), sig_disp));
        dist_factor(counter) = Type(2.0) * dnorm(distances(j), Type(0.0), sig_disp(i));
        fifty_pct(i) = qnorm(Type(0.75), Type(0.0), sig_disp(i)); // want middle of positive half of distn
        pct_at_dist(i) = (pnorm(dist_cutoff, Type(0.0), sig_disp(i)) - Type(0.5)) * Type(2.0);
      } else if(disp_model == 2) { // 2 = exponential
        // dist_factor(counter) = (pexp(distances(j)+site_width, 1/sig_disp) - //F(d+site_width)
        //   pexp(distances(j)-site_width, 1/sig_disp)); //F(d-site_width)
        dist_factor(counter) = dexp(distances(j), pow(sig_disp(i), -1));
        fifty_pct(i) = qexp(Type(0.5), pow(sig_disp(i), -1));
        pct_at_dist(i) = pexp(dist_cutoff, pow(sig_disp(i), -1));
      } else if(disp_model == 3) { // 3 = half-cauchy
        // dist_factor(counter) = (Type(2.0)/PI) * atan((distances(j)+site_width)/sig_disp) - //F(d+site_width)
        //   (Type(2.0)/PI) * atan((distances(j)-site_width)/sig_disp); //F(d-site_width)
        dist_factor(counter) = Type(2.0) / (PI*sig_disp(i) * (1 + pow(distances(j)/sig_disp(i),2)));
        fifty_pct(i) = sig_disp(i)*tan(PI*(Type(0.75) - Type(0.5))); // want middle of positive half of distn
        pct_at_dist(i) = pow(PI, -1) * atan(dist_cutoff/sig_disp(i)) * Type(2.0);
      }
      
      // predicted counts, NLL, outmatrix, etc.
      predcount(counter) = nrel * pow(survival, times(i)) * detectability * dist_factor(counter);
      // obscount(counter) = countmat(j+(i*nsites), k);
      if(count_model == 1) {
        f -= dpois(obscount(counter), predcount(counter), true);
      } else {
        f -= dnbinom2(obscount(counter), predcount(counter), 
                      predcount(counter) + pow(predcount(counter), 2)/overdispersion, true);
      }
      
      outmat(counter, 0) = obscount(counter);
      outmat(counter, 1) = predcount(counter);
      outmat(counter, 2) = times(i);
      outmat(counter, 3) = distances(j);
      outmat(counter, 4) = dist_factor(counter);
      
      counter++;
      // }
    }   
  }
  
  ADREPORT(survival);
  ADREPORT(detectability);
  ADREPORT(exp(log_sig_disp_mu));
  ADREPORT(sig_disp);
  ADREPORT(fifty_pct);
  ADREPORT(pct_at_dist);
  REPORT(outmat);
  return f;
}

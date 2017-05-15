#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_INTEGER(disp_model); // 1 = half-normal, 2 = exponential, 3 = half-cauchy
  // DATA_INTEGER(ncreeks);  // ncreeks is the number of release locations (creeks). Not currently used.
  DATA_INTEGER(nrel);  // nrel is the number of tagged fish released
  DATA_INTEGER(nsites); // the number of trapping sites at each location
  DATA_INTEGER(nperiods); // the number of recap events at each location
  DATA_INTEGER(ntraps); // the number of traps at each site
  DATA_MATRIX(countmat); //(1,nsites*nperiods,1,ntraps) // the number of recap events at each location
  DATA_VECTOR(distances); //(1,nsites) // vector of distances from release site to center of recap site
  DATA_VECTOR(times); //(1,nperiods) // vector of times since release
  // DATA_SCALAR(site_width);

  // PARAMETER(logit_survival);
  // PARAMETER(log_detectability);
  // PARAMETER_VECTOR(log_sig_disp);

  PARAMETER(survival);
  PARAMETER(detectability);
  // PARAMETER_VECTOR(sig_disp);
  PARAMETER(sig_disp);

  // Type survival = exp(logit_survival)/(1+exp(logit_survival));
  // Type detectability = exp(logit_detectability)/(1+exp(logit_detectability));
  // vector<Type> sig_disp(nperiods);
  // for(int period=0; period<nperiods; period++) {
  //   sig_disp(period) = exp(log_sig_disp(period));
  // }

  // Type site_width = full_site_width / Type(2.0);

  vector<Type> predcount(nsites * nperiods * ntraps);
  vector<Type> obscount(nsites * nperiods * ntraps);
  vector<Type> dist_factor(nsites * nperiods * ntraps);
  matrix<Type> outmat(nsites * nperiods * ntraps, 5);
  Type f = 0; // objective function value

  int counter = 0;
  for (int i=0; i<nperiods; i++) {
    for (int j=0; j<nsites; j++) {
      for (int k=0; k<ntraps; k++) {

      if(disp_model == 1) { // 1 = half-normal
        // dist_factor(counter) = Type(2.0) * (pnorm(distances(j)+site_width, Type(0.0), sig_disp) - 
        //     pnorm(distances(j)-site_width, Type(0.0), sig_disp));
        dist_factor(counter) = Type(2.0) * dnorm(distances(j), Type(0.0), sig_disp);
      } else if(disp_model == 2) { // 2 = exponential
        // dist_factor(counter) = (pexp(distances(j)+site_width, 1/sig_disp) - //F(d+site_width)
        //   pexp(distances(j)-site_width, 1/sig_disp)); //F(d-site_width)
        dist_factor(counter) = dexp(distances(j), pow(sig_disp, -1));
      } else if(disp_model == 3) { // 3 = half-cauchy
        // dist_factor(counter) = (Type(2.0)/PI) * atan((distances(j)+site_width)/sig_disp) - //F(d+site_width)
        //   (Type(2.0)/PI) * atan((distances(j)-site_width)/sig_disp); //F(d-site_width)
        dist_factor(counter) = Type(2.0) / (PI*sig_disp * (1 + pow(distances(j)/sig_disp,2)));
      }

     predcount(counter) = nrel * pow(survival, times(i)) * detectability * dist_factor(counter);
     obscount(counter) = countmat(j+(i*nsites), k);
     f -= dpois(obscount(counter), predcount(counter), true);

     outmat(counter, 0) = obscount(counter);
     outmat(counter, 1) = predcount(counter);
     outmat(counter, 2) = times(i);
     outmat(counter, 3) = distances(j);
     outmat(counter, 4) = dist_factor(counter);
     
     counter++;
    }
   }   
  }
 
  ADREPORT(survival);
  ADREPORT(detectability);
  ADREPORT(sig_disp);
  REPORT(outmat);
  return f;
}

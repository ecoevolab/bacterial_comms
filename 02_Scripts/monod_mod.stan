// MONOD EQUATION

functions {
  vector lvfnc(real t,           // time
               vector z,         // state: z[1] = Biomass (OD), z[2] = Substrate
               real r,           // r value (mu_max in the original model)
               real ka,          // affinity constant (k_s in the original model)
               real yield) {     // (Y) 
    
    vector[2] dzdt ;     // vecttor dzdt 
    
    
    // to avoid 0 values 
    real X = fmax(z[1], 1e-6); // Biomass (OD)
    real S = fmax(z[2], 1e-6); // Substrate
    
    dzdt[1] = r * (S / (ka + S)) * X; // microbial growth 
    dzdt[2] = - (1.0 / yield) * dzdt[1];  //substrate consumption 
    
    // note: dzdt[2] = substrate -- is going to decrease by consumption / proportional to generated biomass
    
    }
    return dzdt;
  }
}

data {
  int<lower=0> S;             // Number of REPLICAS 
  int<lower=0> totalobs;       // sum of all the points in the vector
  array[S] int<lower=1> sizes; // how many points there are in each replica (timepoints)
  
  vector[totalobs] y_time;         // Vector list with time measures (plane)
  vector[totalobs] y_obs;         // vector list with population measures (all of them )
  
  real<lower=0> sigma;        // specific value for sigma 
  real rin;                   // for the specific r prior
  real kin;                   // for the specific k prior 
  real kfin;                  // to the to OD value 
  real sdk;                   // standard deviation for k & r respectively 
  real sdr;
}


parameters {
  real<lower=0> r;           // r value 
  real<lower=0> ka;           // ka
  real<lower=kin, upper=kfin> yield;           // yield value 
}

model {
  // Priors
  r ~ lognormal(log(rin), sdr);
  ka ~ lognormal(log(0.05), sdk);
  yield ~ lognormal(log(kin), sdk);

  // Prior for the initial state based on the values
  int point = 1; 
  
  // extract the data for THAT replica 
  for (i in 1:S){                                  // s = number of replicas "do it s times"
   int N_s = sizes[i];                             // to get the size of the n vector (replica 1 = 10 timepoints)
   int N_sim = N_s - 1;

   vector[N_s] time_s = segment(y_time, point, N_s); // create a 10 n vector / time_s / y_time = source, point = initial value of the sample, N_s = size of the vector 
   vector[N_s] obs_s = segment(y_obs, point, N_s);   // extract the specific observations
   
   // extract the initial value 
   vector[2] z0;
   z0[1] = fmax(obs_s[1], 1e-3) ; // z0 size 1 / its going to be obs_s first value (initial sampling value)
   z0[2] = 1.0;   // for normalizing initial substrate at 100% 
   
   real t_start = time_s[1]; // get the first time for the start of the simulation 

   // modify the segment as array 
    array[N_sim] real t_array; // create an array, of size N_sim (N_s - 1 (because we don't want to count the first time))
    
    for (n in 1:N_sim) { // "do it 1-x times"
      t_array[n] = time_s[n + 1];// copy the timeseries to the t_array (for ode_rk45 to function)
      
    }
   
   // solve the ODE with the array 
   array[N_sim] vector[2] z = ode_rk45(lvfnc, z0, t_start, t_array, r, ka, yield);   // solve it 
   
   // likelihood 
   for (n in 1:N_sim) {
    obs_s[n + 1] ~ normal(z[n][1], sigma);
  // fmax compares two numbers and get the higher one - if the ode number is near 0 fmax 
  // replace it with 0.000001 this to get strictly positive numbers 
   }
   
   obs_s[1] ~ normal(z0[1], sigma);

   point = point + N_s;
}
}


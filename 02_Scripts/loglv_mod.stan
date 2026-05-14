
functions {
  vector lvfnc(real t,           // time
               vector z,         // state (OD in this particular case)
               real r,           // r value
               real k) {         // k value 
    
    vector[num_elements(z)] dzdt ;     // vecttor dzdt 
    
    for (j in 1:num_elements(z)) {
      dzdt[j] = r * z[j] * (1 - z[j] / k);
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
}


parameters {
  real<lower=0> r;           // r value 
  real<lower=0> k;           // k value 
}

model {
  // Priors
  r ~ lognormal(0, 0.5);
  k ~ lognormal(0.2, 0.5);

  // Prior for the initial state based on the values
  int point = 1; 
  
  // extract the data for THAT replica 
  for (i in 1:S){                                  // s = number of replicas "do it s times"
   int N_s = sizes[i];                             // to get the size of the n vector (replica 1 = 10 timepoints)
   int N_sim = N_s - 1;

   vector[N_s] time_s = segment(y_time, point, N_s); // create a 10 n vector / time_s / y_time = source, point = initial value of the sample, N_s = size of the vector 
   vector[N_s] obs_s = segment(y_obs, point, N_s);   // extract the specific observations
   
   // extract the initial value 
   vector[1] z0;
   z0[1] = obs_s[1]; // z0 size 1 / its going to be obs_s first value (initial sampling value)
   real t_start = time_s[1]; // get the first time for the start of the simulation 

   // modify the segment as array 
    array[N_sim] real t_array; // create an array, of size N_sim (N_s - 1 (because we don't want to count the first time))
    
    for (n in 1:N_sim) { // "do it 1-x times"
      t_array[n] = time_s[n + 1]; #copy the timeseries to the t_array (for ode_rk45 to function)
    }
   
   // solve the ODE with the array 
array[N_sim] vector[1] z = ode_rk45(lvfnc, z0, t_start, t_array, r, k);   #solve it 
   
   // likelihood 
   for (n in 1:N_sim){
     obs_s[n + 1] ~ normal(z[n][1], sigma);
   }
   
   point = point + N_s; # to get the next vector size and repeat S
  }



}

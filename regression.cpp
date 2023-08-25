
// linear regression
#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() () {
// data:
DATA_VECTOR(x);
DATA_VECTOR(y);

// parameters:
PARAMETER(a); // intercept
PARAMETER(b); // slope
PARAMETER(log_sigma); // log(residual SD)
// we fit sigma on a log scale to keep it > 0

// procedures: (transformed parameters)
Type sigma = exp(log_sigma);

int n = y.size(); // get number of data points to loop over

Type nll = 0.0; // initialize negative log likelihood

for(int i = 0; i < n; i++){ // C++ starts loops at 0!
  // get negative log likelihood (last argument is log = TRUE)
  nll -= dnorm(y[i], a + b * x[i], sigma, true);
}

return nll;
}

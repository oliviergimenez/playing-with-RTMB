// linear regression
#include <TMB.hpp>

//template<class Type>
//matrix<Type> multmat(array<Type> A, matrix<Type> B) {
//  int nrowa = A.rows();
//  int ncola = A.cols(); 
//  int ncolb = B.cols(); 
//  matrix<Type> C(nrowa,ncolb);
//  for (int i = 0; i < nrowa; i++)
//  {
//    for (int j = 0; j < ncolb; j++)
//    {
//      C(i,j) = Type(0);
//      for (int k = 0; k < ncola; k++)
//        C(i,j) += A(i,k)*B(k,j);
//    }
//  }
//  return C;
//}
//
/* implement the vector - matrix product */
template<class Type>
vector<Type> multvecmat(array<Type>  A, matrix<Type>  B) {
  int nrowb = B.rows();
  int ncolb = B.cols(); 
  vector<Type> C(ncolb);
  for (int i = 0; i < ncolb; i++)
  {
    C(i) = Type(0);
    for (int k = 0; k < nrowb; k++){
      C(i) += A(k)*B(k,i);
    }
  }
  return C;
}

template<class Type>
Type objective_function<Type>::operator() () {
  
  // b = parameters
  PARAMETER_VECTOR(b);
  
  // ch = capture-recapture histories (individual format)
  // fc = date of first capture
  // fs = state at first capture
  DATA_IMATRIX(ch);
  DATA_IVECTOR(fc);
  DATA_IVECTOR(fs);
  
  // OBSERVATIONS
  // 0 = non-detected
  // 1 = detected
  //   
  // STATES
  // 1 = alive 
  // 2 = dead
  //   
  // PARAMETERS
  // phi  survival prob.
  // p  detection prob.
  // pi prob. of being in initial state alive
  //   
  // logit link for all parameters

  int km = ch.rows();
  int nh = ch.cols();  
  int npar = b.size();
  vector<Type> par(npar);
  for (int i = 0; i < npar; i++) {
    par(i) = Type(1.0) / (Type(1.0) + exp(-b(i)));
  }
//  Type pi = par(0); // careful, indexing starts at 0 in Rcpp!
  Type phi = par(0);
  Type p = par(1);

  // prob of obs (rows) cond on states (col)
  matrix<Type> B(2,2);
  B(0,0) = Type(1.0)-p;
  B(0,1) = Type(1.0);
  B(1,0) = p;
  B(1,1) = Type(0.0);

  // first encounter
  matrix<Type> BE(2, 2);
  BE(0,0) = Type(0.0);
  BE(0,1) = Type(1.0);
  BE(1,0) = Type(1.0);
  BE(1,1) = Type(0.0);

  // prob of states at t+1 given states at t
  matrix<Type> A(2, 2);
  A(0,0) = phi;
  A(0,1) = Type(1.0)-phi;
  A(1,0) = Type(0.0);
  A(1,1) = Type(1.0);
  
  // init states
  vector<Type> PROP(2);
  PROP(0) = Type(1.0);
  PROP(1) = Type(0.0);
//  PROP(0) = pi;
//  PROP(1) = Type(1.0)-pi;
//  REPORT(PROP);
  
  // likelihood
  Type ll;
  Type nll;
  array<Type> ALPHA(2);
  for (int i = 0; i < nh; i++) {
    int ei = fc(i)-1;
    vector<int> evennt = ch.col(i);
    ALPHA = PROP * vector<Type>(BE.row(fs(i))); // element-wise vector product
    for (int j = ei+1; j < km; j++) {
      ALPHA = multvecmat(ALPHA,A) * vector<Type>(B.row(evennt(j))); // vector matrix product, then element-wise vector product
    }
    ll += log(sum(ALPHA));
  }
  ADREPORT(p);
  ADREPORT(phi);
  nll = -ll;
  return nll;
}


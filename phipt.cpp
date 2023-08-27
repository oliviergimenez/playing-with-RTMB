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
  // p  detection prob. timedep
  // pi prob. of being in initial state alive
  //   
  // logit link for all parameters

  int km = ch.rows();
  int nh = ch.cols();  

  // b = parameters
  PARAMETER(phib); Type phi = Type(1.0) / (Type(1.0) + exp(-phib)); 
  PARAMETER_VECTOR(detb); 
  vector<Type> p(6);
  p(0) = Type(1.0) / (Type(1.0) + exp(-detb(0)));
  p(1) = Type(1.0) / (Type(1.0) + exp(-detb(1)));
  p(2) = Type(1.0) / (Type(1.0) + exp(-detb(2)));
  p(3) = Type(1.0) / (Type(1.0) + exp(-detb(3)));
  p(4) = Type(1.0) / (Type(1.0) + exp(-detb(4)));
  p(5) = Type(1.0) / (Type(1.0) + exp(-detb(5)));
  
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
  Type ll = 0;
  Type nll = 0;
  array<Type> ALPHA(2);
  for (int i = 0; i < nh; i++) {
    int ei = fc(i)-1;
    vector<int> evennt = ch.col(i);
    ALPHA = PROP * vector<Type>(BE.row(fs(i))); // element-wise vector product
    for (int j = ei+1; j < km; j++) {
      
      // prob of obs (rows) cond on states (col)
      matrix<Type> B(2,2);
      B(0,0) = Type(1.0) - p(j-1);
      B(0,1) = Type(1.0);
      B(1,0) = p(j-1);
      B(1,1) = Type(0.0);
      
      ALPHA = multvecmat(ALPHA,A) * vector<Type>(B.row(evennt(j))); // vector matrix product, then element-wise vector product
    }
    ll += log(sum(ALPHA));
  }
  nll = -ll;
  ADREPORT(phi);
  ADREPORT(p);
  return nll;
}


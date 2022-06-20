#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;
using namespace std;

// [[Rcpp::depends(RcppArmadillo)]]

uword getNegative(arma::mat X) {
  uvec q1 = find(X < 0);
  vec B = conv_to<vec>::from(q1);
  return B.n_elem;
}
// [[Rcpp::export]]
double getSum(arma::mat X, arma::mat M) {
  return accu(X) / M.n_rows;
}


arma::vec nnls_col(const mat &A, const subview_col<double> &b, int max_iter = 500, double tol = 1e-6, bool verbose = false)
{
  
  vec mu = -A.t() * b;
  mat H = A.t() * A;
  vec x(A.n_cols), x0(A.n_cols);
  x.fill(0);
  x0.fill(-9999);
  
  int i = 0;
  double tmp;
  while(i < max_iter && max(abs(x - x0)) > tol) 
  {
    x0 = x;
    for (int k = 0; k < A.n_cols; k++) 
    {
      tmp = x[k] - mu[k] / H.at(k,k);
      if (tmp < 0) tmp = 0;
      if (tmp != x[k]) mu += (tmp - x[k]) * H.col(k);
      x[k] = tmp;
    }
    ++i;
  }
  
  return x;
}


arma::mat nnls(mat A, mat b, int max_iter = 500, double tol = 1e-6)
{
  // solving Ax = b, where x and b are both matrices
  if(A.n_rows != b.n_rows)
    throw std::invalid_argument("A and b must have the same number of rows.");
  mat x(A.n_cols, b.n_cols);
  for (int i = 0; i < b.n_cols; i++)
    x.col(i) = nnls_col(A, b.col(i), max_iter, tol);
  
  return x;
}



arma::mat hinge_der_proportions_C__(const arma::mat& H,
                                    const arma::mat& R, double precision_ = 1e-10) {
  int m = H.n_rows;
  int n = H.n_cols;
  
  arma::mat TMP(n, m*m, fill::zeros);
  
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++) {
      if (H(i,j) < 0) {
        for (int k = 1; k < m; k++) {
          TMP(j,k+i*m) = -R(k,j);
        }
      }
    }
  }
  
  return reshape(arma::sum(TMP,0),m,m).t();
  
}

arma::mat hinge_der_basis_C__(const arma::mat& W,
                              const arma::mat& S, 
                              double precision_ = 1e-10) {
  int n = W.n_cols;
  
  arma::mat res(n, n, fill::zeros);
  
  for (int j = 0; j < n; j++) {
    
    arma::vec t = W.col(j); 
    res.col(j) =  arma::sum(-S.cols(find(t < -precision_)),1);
  }
  res.row(0).zeros();
  
  return res;
  
}

double hinge_C__(const arma::mat& X) {
  arma::mat X_(X.n_rows,X.n_cols,fill::zeros);
  double elem_ = 0;
  
  for (int i=0; i<X.n_rows; i++) {
    for (int j=0; j<X.n_cols; j++) {
      elem_ = X(i,j);
      if (elem_<0){
        X_(i,j) = -elem_;
      }
    }
  }
  return accu(X_);
}

List calcErrors(const arma::mat& X,
                const arma::mat& Omega,
                const arma::mat& D_w,
                const arma::mat& D_h,
                const arma::mat& V_row,
                const arma::mat& R,
                const arma::mat& S,
                const double coef_,
                const double coef_der_X,
                const double coef_der_Omega,
                const double coef_hinge_H,
                const double coef_hinge_W,
                const double coef_pos_D_h,
                const double coef_pos_D_w) {
  arma::mat V__ = S * V_row * R.t();
  arma::mat D_w_diag = diagmat(D_w);
  
  double error_ = pow(norm(V__ - Omega * D_w_diag * X,"fro"),2.0);
  double orig_deconv_error = pow(norm(V_row - S.t() * Omega * D_w_diag * X * R,"fro"),2);
  double lambda_error = coef_ * coef_hinge_H * hinge_C__(X * R);
  double beta_error = coef_ * coef_hinge_W * hinge_C__(S.t() * Omega);
  arma::mat A = arma::sum(R,1);
  arma::mat B = arma::sum(S,1);
  double D_h_error = coef_pos_D_h * pow(norm(X.t() * D_h - A,"fro"),2);
  double D_w_error = coef_pos_D_w * pow(norm(Omega * D_w - B,"fro"),2);
  double new_error = error_ + lambda_error + beta_error + D_h_error + D_w_error;
  
  return List::create(Named("error_") = error_,
                      Named("lambda_error") = lambda_error,
                      Named("beta_error") = beta_error,
                      Named("D_h_error") = D_h_error,
                      Named("D_w_error") = D_w_error,
                      Named("new_error") = new_error,
                      Named("orig_deconv_error") = orig_deconv_error);
}

field<mat> derivative_stage1(const arma::mat& X,
                           const arma::mat& Omega,
                           const arma::mat& D_w,
                           const arma::mat& V_row,
                           const arma::mat& R,
                           const arma::mat& S,
                           const arma::vec& splits,
                           const arma::vec& intervals,
                           const double coef_der_X,
                           const double coef_der_Omega,
                           const double coef_hinge_H,
                           const double coef_hinge_W,
                           const double coef_pos_D_h,
                           const double coef_pos_D_w,
                           const int cell_types,
                           const double N,
                           const double M,
                           const int global_iterations,
                           arma::mat& errors_statistics,
                           arma::mat& points_statistics_X,
                           arma::mat& points_statistics_Omega,
                           const int steps = 5) {
  
  arma::mat new_X = X;
  arma::mat new_Omega = Omega;
  arma::mat new_D_w = D_w;
  arma::mat new_D_h = new_D_w * (N/M);
  double coef_;
  
  arma::mat V__ = S * V_row * R.t();
  mat B = join_cols(arma::vectorise(V__),coef_pos_D_w * arma::sum(S,1));
  arma::mat der_X,der_Omega;
  
  for (int itr_= 1; itr_ <= global_iterations; itr_++) {
    int idx_ = intervals(itr_-1);
    
    if (idx_ % 2 == 0) {
      coef_ = splits(idx_-1);
    } else {
      coef_ = 1;
    }
    
    for (int i=0; i < steps; i++) {
      // derivative X
      der_X = -2 * (diagmat(new_D_w) * new_Omega.t() * (V__ - new_Omega * diagmat(new_D_w) * new_X));
      der_X = der_X + coef_ * coef_hinge_H * hinge_der_proportions_C__(new_X * R, R);
      der_X.col(0).zeros();
      
      new_X = new_X - coef_der_X * der_X;
      
      // derivative Omega
      der_Omega = -2 * (V__ - new_Omega * diagmat(new_D_w) * new_X) * new_X.t() * diagmat(new_D_w);
      der_Omega = der_Omega + coef_ * coef_hinge_W * hinge_der_basis_C__(S.t() * new_Omega, S);
      der_Omega.row(0).zeros();
      
      new_Omega = new_Omega - coef_der_Omega * der_Omega;
    }
    
    arma::mat vec_mtx(cell_types*cell_types,cell_types,fill::zeros);
    for (int c=0; c<cell_types; c++) {
      vec_mtx.col(c) = arma::vectorise(new_Omega.col(c) * new_X.row(c));
    }
    arma::mat A = join_cols(vec_mtx, coef_pos_D_w * new_Omega);
    
    
    new_D_w = nnls(A, B);
    new_D_w.elem(find(new_D_w <= 0)).fill(1e-09);
    
    new_D_h = new_D_w * (N/M);
    
    uword neg_props = getNegative(new_X * R);
    uword neg_basis = getNegative(S.t() * new_Omega);
    double sum_ = accu(new_D_w) / V_row.n_rows;
    
    List err_ = calcErrors(new_X,new_Omega,new_D_w, new_D_h,
                           V_row, R, S, coef_, coef_der_X,
                           coef_der_Omega, coef_hinge_H,
                           coef_hinge_W, coef_pos_D_h,
                           coef_pos_D_w);
    errors_statistics.row(itr_) = { err_["error_"],
                                   err_["lambda_error"],
                                   err_["beta_error"],
                                   err_["D_h_error"],
                                   err_["D_w_error"],
                                   err_["new_error"],
                                   err_["orig_deconv_error"],
                                   neg_props,
                                   neg_basis,
                                   sum_};
    points_statistics_X.row(itr_) = new_X.as_row();
    points_statistics_Omega.row(itr_) = new_Omega.as_row();
  }
  
  field<mat> ret_(7,1);
  ret_(0,0) = new_X;
  ret_(1,0) = new_Omega;
  ret_(2,0) = new_D_w;
  ret_(3,0) = new_D_h;
  ret_(4,0) = errors_statistics;
  ret_(5,0) = points_statistics_X;
  ret_(6,0) = points_statistics_Omega;
  
  return ret_;
}

// [[Rcpp::export]]
field<mat> derivative_stage2(const arma::mat& X,
                             const arma::mat& Omega,
                             const arma::mat& D_w,
                             const arma::mat& V_row,
                             const arma::mat& R,
                             const arma::mat& S,
                             const double coef_der_X,
                             const double coef_der_Omega,
                             const double coef_hinge_H,
                             const double coef_hinge_W,
                             const double coef_pos_D_h,
                             const double coef_pos_D_w,
                             const int cell_types,
                             const double N,
                             const double M,
                             const int global_iterations,
                             arma::mat& errors_statistics,
                             const int start_idx,
                             arma::mat& points_statistics_X,
                             arma::mat& points_statistics_Omega) {
  arma::mat new_X = X;
  arma::mat new_Omega = Omega;
  arma::mat new_D_w = D_w;
  arma::mat new_D_h = new_D_w * (N/M);
  
  arma::mat V__ = S * V_row * R.t();
  mat B = join_cols(arma::vectorise(V__),coef_pos_D_w * arma::sum(S,1));
  arma::mat der_X,der_Omega;
  
  for (int itr_= start_idx; itr_ < global_iterations+start_idx; itr_++) {
    
    // derivative X
    der_X = -2 * (diagmat(new_D_w) * new_Omega.t() * (V__ - new_Omega * diagmat(new_D_w) * new_X));
    der_X += coef_hinge_H * hinge_der_proportions_C__(new_X * R, R);
    der_X += coef_pos_D_h * 2 * new_D_h * (new_X.t() * new_D_h - arma::sum(R,1)).t();
    der_X.col(0).zeros();
    
    new_X = new_X - coef_der_X * der_X;
    
    arma::mat vec_mtx(cell_types*cell_types,cell_types,fill::zeros);
    for (int c=0; c<cell_types; c++) {
      vec_mtx.col(c) = arma::vectorise(new_Omega.col(c) * new_X.row(c));
    }
    arma::mat A = join_cols(vec_mtx, coef_pos_D_w * new_Omega);
    
    
    new_D_w = nnls(A, B);
    new_D_w.elem(find(new_D_w <= 0)).fill(1e-09);
    
    new_D_h = new_D_w * (N/M);
    
    // derivative Omega
    der_Omega = -2 * (V__ - new_Omega * diagmat(new_D_w) * new_X) * new_X.t() * diagmat(new_D_w);
    der_Omega += coef_hinge_W * hinge_der_basis_C__(S.t() * new_Omega, S);
    der_Omega += coef_pos_D_w * 2 * (new_Omega*new_D_w-arma::sum(S,1)) * D_w.t();
    der_Omega.row(0).zeros();
    
    new_Omega = new_Omega - coef_der_Omega * der_Omega;
    
    
    vec_mtx.fill(fill::zeros);
    A.fill(fill::zeros);
    
    for (int c=0; c<cell_types; c++) {
      vec_mtx.col(c) = arma::vectorise(new_Omega.col(c) * new_X.row(c));
    }
    A = join_cols(vec_mtx, coef_pos_D_w * new_Omega);
    
    
    new_D_w = nnls(A, B);
    new_D_w.elem(find(new_D_w <= 0)).fill(1e-09);
    
    new_D_h = new_D_w * (N/M);
    
    uword neg_props = getNegative(new_X * R);
    uword neg_basis = getNegative(S.t() * new_Omega);
    double sum_ = accu(new_D_w) / V_row.n_rows;
    
    List err_ = calcErrors(new_X,new_Omega,new_D_w, new_D_h,
                           V_row, R, S, 1, coef_der_X,
                           coef_der_Omega, coef_hinge_H,
                           coef_hinge_W, coef_pos_D_h,
                           coef_pos_D_w);
    errors_statistics.row(itr_) = { err_["error_"],
                                    err_["lambda_error"],
                                        err_["beta_error"],
                                            err_["D_h_error"],
                                                err_["D_w_error"],
                                                    err_["new_error"],
                                                        err_["orig_deconv_error"],
                                                            neg_props,
                                                            neg_basis,
                                                            sum_};
    points_statistics_X.row(itr_) = new_X.as_row();
    points_statistics_Omega.row(itr_) = new_Omega.as_row();
    
    
  }
  
  field<mat> ret_(7,1);
  ret_(0,0) = new_X;
  ret_(1,0) = new_Omega;
  ret_(2,0) = new_D_w;
  ret_(3,0) = new_D_h;
  ret_(4,0) = errors_statistics;
  ret_(5,0) = points_statistics_X;
  ret_(6,0) = points_statistics_Omega;
  
  return ret_;
  
}

// [[Rcpp::export]]
List run_optimization( const arma::mat& X,
                       const arma::mat& Omega,
                       const arma::mat& D_w,
                       const arma::mat& V_row,
                       const arma::mat& R,
                       const arma::mat& S,
                       const arma::vec& splits,
                       const arma::vec& intervals,
                       const double coef_der_X,
                       const double coef_der_Omega,
                       const double coef_hinge_H,
                       const double coef_hinge_W,
                       const double coef_pos_D_h,
                       const double coef_pos_D_w,
                       const int cell_types,
                       const double N,
                       const double M,
                       const int global_iterations,
                       const bool startWithInit ) {
  arma::mat errors_statistics;
  arma::mat points_statistics_X;
  arma::mat points_statistics_Omega;
  
  int start_idx, rows_iterations;
  mat new_X = X;
  mat new_Omega = Omega;
  mat new_D_w = D_w;
  mat new_D_h = D_w * (N/M);
  
  if (startWithInit) {
    rows_iterations = 2 * global_iterations + 1;
  } else {
    rows_iterations = global_iterations + 1;
  }
  errors_statistics.set_size(rows_iterations, 10);
  errors_statistics.fill(fill::zeros);
  
  points_statistics_X.set_size(rows_iterations, cell_types*cell_types);
  points_statistics_X.fill(fill::zeros);
  
  points_statistics_Omega.set_size(rows_iterations, cell_types*cell_types);
  points_statistics_Omega.fill(fill::zeros);
  
  List err_ = calcErrors(new_X, new_Omega, new_D_w, new_D_h,
                         V_row, R, S, 1, coef_der_X,
                         coef_der_Omega, coef_hinge_H,
                         coef_hinge_W, coef_pos_D_h,
                         coef_pos_D_w);
  
  uword neg_props = getNegative(new_X * R);
  uword neg_basis = getNegative(S.t() * new_Omega);
  double sum_ = accu(new_D_w) / V_row.n_rows;
  
  errors_statistics.row(0) = { err_["error_"],
                               err_["lambda_error"],
                               err_["beta_error"],
                               err_["D_h_error"],
                               err_["D_w_error"],
                               err_["new_error"],
                               err_["orig_deconv_error"],
                                   neg_props,
                                   neg_basis,
                                   sum_};
  
  points_statistics_X.row(0) = new_X.as_row();
  points_statistics_Omega.row(0) = new_Omega.as_row();
  
  if (startWithInit) {
    field<mat> stage1 = derivative_stage1(new_X,new_Omega,new_D_w,V_row,
                             R,S,splits,intervals,
                             coef_der_X,coef_der_Omega,
                             coef_hinge_H,coef_hinge_W,
                             coef_pos_D_h,coef_pos_D_w,
                             cell_types,N,M,global_iterations,
                             errors_statistics,
                             points_statistics_X,
                             points_statistics_Omega);
    start_idx = global_iterations+1;
    new_X = stage1(0,0);
    new_Omega = stage1(1,0);
    new_D_w = stage1(2,0);
    new_D_h = stage1(3,0);
    
  } else {
    start_idx = 1;
  }
  
  field<mat> stage2 = derivative_stage2(new_X,new_Omega,new_D_w,V_row,
                                        R,S,coef_der_X,coef_der_Omega,
                                        coef_hinge_H,coef_hinge_W,
                                        coef_pos_D_h,coef_pos_D_w,
                                        cell_types, N, M, global_iterations,
                                        errors_statistics, start_idx,
                                        points_statistics_X,
                                        points_statistics_Omega);
  
  new_X = stage2(0,0);
  new_Omega = stage2(1,0);
  new_D_w = stage2(2,0);
  new_D_h = stage2(3,0);
  
  return List::create(Named("new_X") = new_X,
                      Named("new_Omega") = new_Omega,
                      Named("new_D_w") = new_D_w,
                      Named("new_D_h") = new_D_h,
                      Named("errors") = errors_statistics,
                      Named("points_X") = points_statistics_X,
                      Named("points_Omega") = points_statistics_Omega);
  
  
}





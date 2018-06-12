#include <RcppArmadillo.h>
#include <Rmath.h>
//[[Rcpp::depends("RcppArmadillo")]]


using namespace Rcpp;
using namespace arma;

//[[Rcpp::export]]
arma::vec colSumsCpp(const arma::mat & X){
  int nCols = X.n_cols;
  arma::vec out(nCols);
  for(int i = 0; i < nCols; i++){
    out(i) = sum(X.col(i));
  }
  return(out);
}

//[[Rcpp::export]]
arma::vec rowSumsCpp(const arma::mat & X){
  int nRows = X.n_rows;
  arma::vec out(nRows);
  for(int i = 0; i < nRows; i++){
    out(i) = sum(X.row(i));
  }
  return(out);
}

//[[Rcpp::export]]
arma::vec standCpp (arma::vec x) {
  double mean_x = 1.0/x.n_elem * sum (x);
  double s2_x = 1.0/x.n_elem * sum (square (x - mean_x));
  arma::vec stand_x = (x - mean_x)/sqrt (s2_x);
  return (stand_x);
}

//[[Rcpp::export]]
arma::vec wstandCpp(arma::vec x, arma::vec w) {
  w = w/sum (w);
  arma::mat W = diagmat (w);
  arma::colvec I = ones<colvec>(x.n_elem, 1);
  double wmean_x = as_scalar (I.t() * W * x);  // the operation returns matrix, eventhought it should be a double
  double ws2_x = as_scalar ((x-wmean_x*I).t() * W * (x - wmean_x * I));
  arma::vec wstand_x = (x - wmean_x)/sqrt (ws2_x); 
  return (wstand_x);
}

//[[Rcpp::export]]
double wmeanCpp (arma::vec x, arma::vec w) {
  w = w/sum (w);
  arma::mat W = diagmat (w);
  arma::colvec I = ones<colvec>(x.n_elem, 1);
  double wmean_x = as_scalar (I.t() * W * x);  // the operation returns matrix, eventhought it should be a double
  return (wmean_x);
}

//[[Rcpp::export]]
double wsdCpp (arma::vec x, arma::vec w){
  w = w/sum (w);
  arma::mat W = diagmat (w);
  arma::colvec I = ones<colvec>(x.n_elem, 1);
  double wmean_x = as_scalar (I.t() * W * x);
  double wsd2_x = as_scalar ((x-wmean_x*I).t() * W * (x - wmean_x * I));
  double wsd_x = pow (wsd2_x, 0.5);
  return (wsd_x);
}  

//[[Rcpp::export]]
arma::mat wcenterCpp (arma::mat x, arma::vec w){
  //w = w/sum (w);
  int nr = x.n_rows, nc = x.n_cols;  
  if (w.n_elem != (unsigned)nr)
    throw std::string("weights 'w' and data do not match");
  arma::mat x_wc (nr, nc);
  for (int co = 0; co < nc; co++){
    arma::vec x_co = x.col(co);
    double w_mean = wmeanCpp (x_co, w);
    x_wc.col (co) = (x_co - w_mean) % sqrt (w);  // it took me ages to figure out that I should multiply the result by sqrt of w! This is not standard implementation of weighted centering
    }
  return (x_wc);
}

//[[Rcpp::export]]
double corCpp (arma::vec x, arma::vec y) {
  arma::vec x_stand = standCpp (x);
  arma::vec y_stand = standCpp (y);
  double r = 1.0/x.n_elem * sum (x_stand % y_stand);
  return (r);
}

//[[Rcpp::export]]
double wcorCpp (arma::vec x, arma::vec y, arma::vec w) {
  w = w / sum (w);
  arma::mat W = diagmat (w);
  arma::vec x_wstand = wstandCpp (x, w);
  arma::vec y_wstand = wstandCpp (y, w);
  double r_w = as_scalar (x_wstand.t() * W * y_wstand);
  return (r_w);
}

//[[Rcpp::export]]
arma::mat cwmCpp (arma::mat L, arma::mat t, bool wstand) {
  arma::mat P = L/accu (L);
  arma::mat c(L.n_rows, t.n_cols);
  
  int t_ncol = t.n_cols;
  if (wstand){ // standardized version
    for (int co = 0; co < t_ncol; co++){
      arma::vec t_co = t.col (co);
      arma::vec t_co_temp = t_co.elem (find_finite (t_co));
      arma::mat P_temp = P.cols (find_finite (t_co));
      arma::vec w_s = colSumsCpp (P_temp);
      arma::vec t_co_temp_ws = wstandCpp (t_co_temp, w_s);
      arma::vec w_n = rowSumsCpp (P_temp);
      if (any (w_n == 0)) {
        arma::mat W_n = diagmat (w_n.elem (find (w_n > 0)));
        arma::mat P_temp2 = P_temp.rows (find (w_n > 0));
        arma::vec c_col (P.n_rows);
        c_col.fill (NA_REAL);
        c_col.elem (find (w_n > 0)) = inv (W_n) * P_temp2 * t_co_temp_ws;
        c.col (co) = c_col;
      } else {
        arma::mat W_n = diagmat (w_n);
        c.col (co) = inv (W_n) * P_temp * t_co_temp_ws;
        }
      }
  } else { //non-standardized version
      for (int co = 0; co < t_ncol; co++){
        arma::vec t_co = t.col (co);
        arma::vec t_co_temp = t_co.elem (find_finite (t_co));
        arma::mat P_temp = P.cols (find_finite (t_co));
        arma::vec w_n = rowSumsCpp (P_temp);
        if (any (w_n == 0)) {
          arma::mat W_n = diagmat (w_n.elem (find (w_n > 0)));
          arma::mat P_temp2 = P_temp.rows (find (w_n > 0));
          arma::vec c_col (P.n_rows);
          c_col.fill (NA_REAL);
          c_col.elem (find (w_n > 0)) = inv (W_n) * P_temp2 * t_co_temp;
          c.col (co) = c_col;
        } else {
          arma::mat W_n = diagmat (w_n);
          c.col (co) = inv (W_n) * P_temp * t_co_temp;
          }
        }
     }
  return (c);
} 

//[[Rcpp::export]]
arma::mat sncCpp (arma::mat L, arma::mat e, bool wstand) {
  arma::mat P = L/accu (L);
  arma::mat u(L.n_cols, e.n_cols);
  
  int e_ncol = e.n_cols;
  if (wstand){ // standardized version
    for (int co = 0; co < e_ncol; co++){
      arma::vec e_co = e.col (co);
      arma::vec e_co_temp = e_co.elem (find_finite (e_co));
      arma::mat P_temp = P.cols (find_finite (e_co));
      arma::vec w_s = colSumsCpp (P_temp);
      arma::vec w_n = rowSumsCpp (P_temp);
      arma::vec e_co_temp_ws = wstandCpp (e_co_temp, w_n);
      
      if (any (w_s == 0)) {
        arma::mat W_s = diagmat (w_s.elem (find (w_s > 0)));
        arma::mat P_temp2 = P_temp.cols (find (w_s > 0));
        arma::vec u_col (P.n_cols);
        u_col.fill (NA_REAL);
        u_col.elem (find (w_s > 0)) = inv (W_s) * P_temp2.t() * e_co_temp_ws;
        u.col (co) = u_col;
      } else {
        arma::mat W_s = diagmat (w_s);
        u.col (co) = inv (W_s) * P_temp.t() * e_co_temp_ws;
      }
    }
  } else { //non-standardized version
    for (int co = 0; co < e_ncol; co++){
      arma::vec e_co = e.col (co);
      arma::vec e_co_temp = e_co.elem (find_finite (e_co));
      arma::mat P_temp = P.rows (find_finite (e_co));
      arma::vec w_s = colSumsCpp (P_temp);
      if (any (w_s == 0)) {
        arma::mat W_s = diagmat (w_s.elem (find (w_s > 0)));
        arma::mat P_temp2 = P_temp.cols (find (w_s > 0));
        arma::vec u_col (P.n_cols);
        u_col.fill (NA_REAL);
        u_col.elem (find (w_s > 0)) = inv (W_s) * P_temp2.t() * e_co_temp;
        u.col (co) = u_col;
      } else {
        arma::mat W_s = diagmat (w_s);
        u.col (co) = inv (W_s) * P_temp.t() * e_co_temp;
      }
    }
  }
  return (u);
} 

//[[Rcpp::export]]
arma::vec sncCpp2 (arma::mat L, arma::vec e) {
  arma::mat P = L/accu (L);
  arma::vec w_s = colSumsCpp (P);
  arma::mat W_s = diagmat (w_s);
  arma::vec u = inv (W_s) * P.t() * e;
  return (u);
}

//[[Rcpp::export]]
double ca_eig1_sqrt (arma::mat Y){
  arma::vec fi_ = rowSumsCpp (Y); 
  arma::vec f_j = colSumsCpp (Y);
  double f_ = sum (fi_);
  arma::vec pi_ = fi_/f_;
  arma::vec p_j = f_j/f_;
  arma::mat E = (fi_ * f_j.t())/f_;
  arma::mat Qbar = (Y - E) % pow (E, -0.5) / sqrt (f_);
  arma::vec svd_res = svd (Qbar);
  double eig1 = svd_res (0);
  return (eig1);
}

//[[Rcpp::export]]
double fourth (arma::vec e, arma::mat L, arma::vec t) {
  arma::mat P = L/accu (L);
  arma::vec w_n = rowSumsCpp (P);
  arma::vec w_s = colSumsCpp (P);
  arma::vec e_wn = wstandCpp (e, w_n);
  arma::vec t_ws = wstandCpp (t, w_s);
  double r_f = as_scalar (e_wn.t() * P * t_ws);
  return (r_f);
}

//[[Rcpp::export]]
arma::mat stand_tot (arma::mat L) {
  for (int i = 0; (unsigned)i < L.n_rows; i++){
    L.row(i) /= sum (L.row (i));
  }
  return (L);
}

//[[Rcpp::export]]
int count_if(LogicalVector x) {
  int counter = 0;
  for(int i = 0; i < x.size(); i++) {
    if(x[i] == TRUE) {
      counter++;
    }
  }
  return counter;
}

//[[Rcpp::export]]
arma::mat submat_nonzero(arma::mat X, arma::vec T) {
  arma::mat Xmat(X.begin(), X.n_rows, X.n_cols, false);
  arma::colvec tIdx(T.begin(), T.size(), false); 
  arma::mat y = Xmat.rows(find(tIdx > 0));
  return y;
}


//[[Rcpp::export]]
bool is_in (CharacterVector x, CharacterVector table){
  bool is_in_temp = FALSE;
  for (int i = 0; i < x.size (); i++){
    for (int j = 0; j < table.size (); j++){
      if (x(i) == table(j)) is_in_temp = TRUE;
      };
  };
  return is_in_temp;
}

//[[Rcpp::export]]
arma::uvec complete_cases (arma::mat X){
  int nrows = X.n_rows;
  arma::uvec keep(nrows, fill::zeros);
  for (int row = 0; row < nrows; row++) 
    if (is_finite(X.row(row))) keep(row) = 1;
    return (keep);
}

//[[Rcpp::export]] 
List rm_missing_eLt (arma::mat e, arma::mat L, arma::mat t) {
  arma::uvec e_complcases = complete_cases (e);
  arma::uvec t_complcases = complete_cases (t);
  arma::uvec c_hasspecies = rowSumsCpp(L.cols (find (t_complcases == 1))) > 0;  // finds rows which have species after removing missing t
  arma::mat L_temp = L.submat (find (e_complcases == 1 && c_hasspecies == 1), find (t_complcases == 1));
  arma::mat t_temp = t.rows (find (t_complcases == 1));
  arma::mat e_temp = e.rows (find (e_complcases == 1 && c_hasspecies == 1));
  return List::create (
      _["e"] = e_temp,
      _["L"] = L_temp,
      _["t"] = t_temp,
      _["n_sit"] = L_temp.n_rows,
      _["n_spe"] = L_temp.n_cols,
      _["e_missing"] = find (e_complcases == 0) + 1,
      _["t_missing"] = find (t_complcases == 0) + 1,
      _["c_missing"] = find (c_hasspecies == 0) + 1
  );
}

//[[Rcpp::export]]
List test_cwm_cor (arma::vec e, arma::mat L, arma::vec t, CharacterVector test, int perm, bool wcor, bool wstand) {
  
  List eLt = rm_missing_eLt (e, L, t);  // removing rows and columns in L with missing e and/or t, resp; removing rows with not enough species to calculate
  arma::mat eLt_e = eLt["e"];
  arma::vec e_temp = arma::vectorise (eLt_e);
  arma::mat L_temp = eLt["L"];
  arma::mat eLt_t = eLt["t"];
  arma::vec t_temp = arma::vectorise (eLt_t);
  arma::vec w_n;
  arma::vec w_s;
  if (wcor || wstand) w_n = rowSumsCpp (L_temp);  // rowsums as weights will be used in case wcor == TRUE (w in cor) or wstand == TRUE (w for t)
  if (wstand) w_s = colSumsCpp (L_temp);  // colsums as weights will be used in case wstand == TRUE (w for e)
  arma::vec cwm = cwmCpp (L_temp, t_temp, wstand);
  double no_samples = L_temp.n_rows;
  double r_obs;
  if (wstand) e_temp = wstandCpp (e_temp, w_n);
  if (wcor) r_obs = wcorCpp (cwm, e_temp, w_n); else r_obs = corCpp (cwm, e_temp);
  double t_obs = r_obs*sqrt ((e_temp.size () - 2)/(1 - pow (r_obs, 2.0)));  // calculates t-value according to Student's t formula
  double P_par = R::pt (t_obs, no_samples-2, TRUE, FALSE);
  P_par = 2*(min (NumericVector::create (P_par, 1-P_par)));
  double P_row = NA_REAL, P_col = NA_REAL, P_max = NA_REAL;

  if (is_in (CharacterVector::create ("standard", "rowbased", "max"), test)){
    double r_exp_sta;
    arma::vec t_exp_sta (perm+1);
    for (int nperm = 0; nperm < perm; nperm++){
      arma::mat e_rand = shuffle (e_temp);
      if (wstand) e_rand = wstandCpp (e_rand, w_n);
      if (wcor){r_exp_sta = wcorCpp (cwm, e_rand, w_n);} else {r_exp_sta = corCpp (cwm, e_rand);}
      t_exp_sta (nperm) = r_exp_sta * sqrt ((e_rand.size () - 2)/(1 - pow (r_exp_sta, 2.0)));
    }
    t_exp_sta (perm) = t_obs; //observed value is put on the last place of the vector
    P_row = sum (abs (t_exp_sta) >= std::abs (t_obs))/(perm + 1.0);
  };

  if (is_in (CharacterVector::create ("modified", "colbased", "max"), test)){
    double r_exp_mod;
    arma::vec t_exp_mod (perm+1);
    for (int nperm = 0; nperm < perm; nperm++){
      arma::vec cwm_rand = cwmCpp (L_temp, shuffle (t_temp), wstand);
      if (wcor) {r_exp_mod = wcorCpp (cwm_rand, e_temp, w_n);} else {r_exp_mod = corCpp (cwm_rand, e_temp);}
      t_exp_mod (nperm) = r_exp_mod*sqrt ((e_temp.size () - 2)/(1 - pow (r_exp_mod, 2.0)));
    }
    t_exp_mod (perm) = t_obs; //observed values is put on the last place of the vector
    P_col = sum (abs (t_exp_mod) >= std::abs (t_obs))/(perm + 1.0);
  };

  if (is_in ("max", test)){
    P_max = max (NumericVector::create(P_col, P_row));
  };

  return List::create (
    _["r"] = r_obs,
    _["t"] = t_obs,
    _["n_sit"] = eLt["n_sit"],
    _["n_spe"] = eLt["n_spe"],
    _["P_par"] = P_par,
    _["P_row"] = P_row,
    _["P_col"] = P_col,
    _["P_max"] = P_max,
    _["e_missing"] = eLt["e_missing"],
    _["t_missing"] = eLt["t_missing"],
    _["c_missing"] = eLt["c_missing"]
  );
}

//[[Rcpp::export]]
List test_fourthCpp (arma::vec e, arma::mat L, arma::vec t, CharacterVector test, int perm) {
  
  List eLt = rm_missing_eLt (e, L, t);  // removing rows and columns in L with missing e and/or t, resp; removing rows with not enough species to calculate
  arma::mat eLt_e = eLt["e"];
  arma::vec e_temp = arma::vectorise (eLt_e);
  arma::mat L_temp = eLt["L"];
  arma::mat eLt_t = eLt["t"];
  arma::vec t_temp = arma::vectorise (eLt_t);
  arma::vec w_n;
  double r_f_obs = fourth (e_temp, L_temp, t_temp);
  double r_fc_obs = r_f_obs/ca_eig1_sqrt(L_temp);
  double P_row = NA_REAL, P_col = NA_REAL, P_max = NA_REAL;
  arma::mat e_rand;
  
  if (is_in (CharacterVector::create ("rowbased", "max"), test)){
    arma::vec r_f_exp_sta (perm+1);
    for (int nperm = 0; nperm < perm; nperm++){
      r_f_exp_sta (nperm) = fourth (shuffle (e_temp), L_temp, t_temp);
    }
    r_f_exp_sta (perm) = r_f_obs; //observed value is put on the last place of the vector
    P_row = sum (abs (r_f_exp_sta) >= std::abs (r_f_obs))/(perm + 1.0);
  };
  
  if (is_in (CharacterVector::create ("colbased", "max"), test)){
    arma::vec r_f_exp_mod (perm+1);
    for (int nperm = 0; nperm < perm; nperm++){
      r_f_exp_mod (nperm) = fourth (e_temp, L_temp, shuffle (t_temp));
    }
    r_f_exp_mod (perm) = r_f_obs; //observed values is put on the last place of the vector
    P_col = sum (abs (r_f_exp_mod) >= std::abs (r_f_obs))/(perm + 1.0);
  };
  
  if (is_in ("max", test)){
    P_max = max (NumericVector::create(P_col, P_row));
  };
  
  return List::create (
      _["r_fc"] = r_f_obs,
      _["r_ch"] = r_fc_obs,
      _["n_sit"] = eLt["n_sit"],
      _["n_spe"] = eLt["n_spe"],
      _["P_row"] = P_row,
      _["P_col"] = P_col,
      _["P_max"] = P_max,
      _["e_missing"] = eLt["e_missing"],
      _["t_missing"] = eLt["t_missing"],
      _["c_missing"] = eLt["c_missing"]
  );
}


// Extended version of fastLm from RcppArmadillo, calculates also F-value and r2
//[[Rcpp::export]]
List fastLm_cwm (const arma::mat& X, const arma::colvec& y, bool stderr_incl) {
  int n = X.n_rows, k = X.n_cols;
  
  arma::colvec coef = solve(X, y);    // fit model y ~ X
  arma::colvec res  = y - X*coef;           // residuals
  
  // std.errors of coefficients
  arma::colvec std_err;
  if (stderr_incl){
    double s2 = std::inner_product(res.begin(), res.end(), res.begin(), 0.0)/(n - k);
    std_err = arma::sqrt(s2 * arma::diagvec(arma::pinv(arma::trans(X)*X))); 
  }
  
  
  arma::colvec fit = X*coef;                // fitted values
  double mss = sum (pow (fit - mean (fit), 2.0));    // sum of squares   
  double rss = sum (pow (res - 0, 2.0));              // residual sum of squares  NEED TO SOLVE THE PROBLEMS WITH DOUBLE vs MAT conversion here!
  double rsq = mss/(mss + rss);           // r2
  int mdf = k - 1;
  int rdf = n - k;      // residual number of degrees of freedom, should be n - no.vars - 1 (but X contains 1 column with 1's for intercept
  double F = (mss/mdf)/(rss/rdf);
  return List::create(_["coefficients"] = coef,
                      _["stderr"] = std_err,
                      _["df.residual"]  = rdf,
                      _["df.model"] = mdf,
                      _["F.value"] = F,
                      _["R.squared"] = rsq
                      );
}


//[[Rcpp::export]]
List test_cwm_lm (arma::mat e, arma::mat L, arma::mat t, CharacterVector test, const char* dependence, int perm, bool wstand) {
  
  List eLt = rm_missing_eLt (e, L, t);  // removing rows and columns in L with missing e and/or t, resp; removing rows with not enough species to calculate
  arma::mat e_temp = eLt["e"];
  arma::mat L_temp = eLt["L"];
  arma::mat t_temp = eLt["t"];
  int n = L_temp.n_rows, p;

  arma::mat cwm = cwmCpp (L_temp, t_temp, wstand);

  List lm_obs;
  if (!strncmp (dependence, "cwm ~ env", 1)) {
   lm_obs = fastLm_cwm (join_rows (ones (e_temp.n_rows), e_temp), cwm, TRUE);
   p = e_temp.n_cols;
   }
   else {
    lm_obs = fastLm_cwm (join_rows (ones (cwm.n_rows), cwm), e_temp, TRUE); 
    p = t_temp.n_cols;
  }

  double rsq_obs = lm_obs["R.squared"];
  double rsq_adj = 1.0 - ((n - 1.0)/(n - p - 1.0))*(1.0 - rsq_obs);  // Calculate rsq_adj using Ezekiel's formula
  double rsq_adj_mod = NA_REAL; //rsq_adj_sta = NA_REAL
  double F_obs = lm_obs["F.value"];
  double P_par = R::pf (F_obs, lm_obs["df.model"], lm_obs["df.residual"], FALSE, FALSE);
  double P_row = NA_REAL, P_col = NA_REAL, P_max = NA_REAL;
 
  if (is_in (CharacterVector::create ("standard", "rowbased", "max"), test)){
    arma::vec rsq_exp_sta (perm);
    arma::vec F_exp_sta (perm+1);
    for (int nperm = 0; nperm < perm; nperm++){
      arma::mat e_rand = shuffle (e_temp, 0);  // rows of e_temp are shuffled
      List lm_exp;
      if (!strncmp (dependence, "cwm ~ env", 1)) {
        lm_exp = fastLm_cwm (join_rows (ones (e_rand.n_rows), e_rand), cwm, FALSE);
      } 
      else {
        lm_exp = fastLm_cwm (join_rows (ones (cwm.n_rows), cwm), e_rand, FALSE);
      }
      rsq_exp_sta (nperm) = lm_exp["R.squared"];  // r.sq for calculation of permutation-based adjusted r.sq
      F_exp_sta (nperm) = lm_exp["F.value"];
    }
    F_exp_sta (perm) = F_obs; //observed values is put on the last place of the vector
    P_row = sum (abs (F_exp_sta) >= std::abs (F_obs))/(perm + 1.0);
    //rsq_adj_sta = 1.0 - (1.0/(1.0-mean (rsq_exp_sta))*(1.0-rsq_obs));  // Calculate rsq_adj_sta using values from standard permutation test
  };
  
  if (is_in (CharacterVector::create ("modified", "colbased", "max"), test)){
    arma::vec rsq_exp_mod (perm);  // for adjusted R2
    arma::vec F_exp_mod (perm+1);
    for (int nperm = 0; nperm < perm; nperm++){
      arma::mat cwm_rand = cwmCpp (L_temp, shuffle (t_temp, 0), wstand);
      List lm_exp;
      if (!strncmp (dependence, "cwm ~ env", 1)) {  
        lm_exp = fastLm_cwm (join_rows (ones (e_temp.n_rows), e_temp), cwm_rand, FALSE);
      }
      else {
        lm_exp = fastLm_cwm (join_rows (ones (cwm_rand.n_rows), cwm_rand), e_temp, FALSE);
      }
      rsq_exp_mod [nperm] = lm_exp["R.squared"];  // r.sq for calculation of permutation-based r.sq
      F_exp_mod [nperm] = lm_exp["F.value"];
    }
    F_exp_mod (perm) = F_obs; //observed values is put on the last place of the vector
    P_col = sum (abs (F_exp_mod) >= std::abs (F_obs))/(perm + 1.0);
    rsq_adj_mod = 1.0 - (1.0/(1.0-mean (rsq_exp_mod))*(1.0-rsq_obs));
  };
  
  if (is_in ("max", test)){
    P_max = max (NumericVector::create(P_col, P_row));
  };
  
  return List::create (
      _["coef"] = lm_obs["coefficients"],
      _["stderr"] = lm_obs["stderr"],                   
      _["r2"] = rsq_obs,
      _["r2adj"] = rsq_adj,
      //_["r2adj_sta"] = rsq_adj_sta,  # R2 adjusted by standard perm test
      _["r2adj_mod"] = rsq_adj_mod,
      _["n_sit"] = eLt["n_sit"],
      _["n_spe"] = eLt["n_spe"],
      _["F"] = F_obs,
      _["P_par"] = P_par,
      _["P_row"] = P_row,
      _["P_col"] = P_col,
      _["P_max"] = P_max,
      _["e_missing"] = eLt["e_missing"],
      _["t_missing"] = eLt["t_missing"],
      _["c_missing"] = eLt["c_missing"]
  
);
}

//[[Rcpp::export]]
List lm_vectorfitCpp (arma::mat X, arma::mat P){
  arma::mat coef = arma::solve (X, P);
  arma::mat P_fit = X * coef;
  arma::mat pow_cor = pow (cor (P_fit, P), 2.0);
  arma::vec r = pow_cor.diag ();
  r.replace (NA_REAL, 0);
  return List::create (
      _["coef"] = coef,
      _["r"] = r);
}

//[[Rcpp::export]]
List vectorfit_cwmCpp (arma::mat X, arma::mat L, arma::mat t, bool wstand, arma::colvec w, int perm, CharacterVector test){
  arma::mat C = cwmCpp (L, t, wstand);
  arma::mat Xw = wcenterCpp (X, w);
  arma::mat Cw = wcenterCpp (C, w);
  //int nc = X.n_cols;
  int np = C.n_cols;
  List lmv = lm_vectorfitCpp (Xw, Cw);
  arma::vec r_obs = lmv["r"], P_row (np), P_col (np), P_max (np), P (np);

  if (is_in (CharacterVector::create ("standard", "rowbased", "max"), test) && perm > 0){
    arma::mat r_row (np, perm);
    for (int p = 0; p < perm; p++){
      arma::mat C_r = shuffle (C);
      arma::mat Cw_r = wcenterCpp (C_r, w);
      arma::vec lmv_r = lm_vectorfitCpp(Xw, Cw_r)["r"];
      r_row.col (p) = lmv_r;
    }
    for (int p = 0; p < perm; p++){
      for (int ro = 0; ro < np; ro++){
        if (r_row (ro, p) >= r_obs(ro)) r_row (ro, p) = 1.0; else r_row (ro, p) = 0.0;
      }
    }
    P_row = (rowSumsCpp (r_row) + 1)/(perm+1);
  }
 
 if (is_in (CharacterVector::create ("modified", "colbased", "max"), test) && perm > 0){
   arma::mat r_col (np, perm);
   for (int p = 0; p < perm; p++){
     //if (Progress::check_abort()) return -1.0; // http://gallery.rcpp.org/articles/using-rcppprogress/ 
     arma::mat C_r = cwmCpp (L, shuffle (t), wstand);
     arma::mat Cw_r = wcenterCpp (C_r, w);
     arma::vec lmv_r = lm_vectorfitCpp(Xw, Cw_r)["r"];
     r_col.col (p) = lmv_r;
   }
   for (int p = 0; p < perm; p++){
     for (int ro = 0; ro < np; ro++){
       if (r_col (ro, p) >= r_obs(ro)) r_col (ro, p) = 1.0; else r_col (ro, p) = 0.0;
     }
   }
   P_col = (rowSumsCpp (r_col) + 1)/(perm+1);
 }
 
  if (is_in ("max", test) && perm > 0){
    P_max = max (P_row, P_col);
  }
  
  if ((perm > 0) && (is_in ("max", test))) P = P_max; else
    if ((perm > 0) && (is_in (CharacterVector::create ("modified", "colbased"), test))) P = P_col; else
      if ((perm > 0) && (is_in (CharacterVector::create ("standard", "rowbased"), test))) P = P_row;
  
  arma::mat lmv_coef = lmv["coef"];  
  return List::create (
      _["coef"] = lmv_coef,
      _["heads"] = normalise (lmv_coef, 2),
      _["r"] = lmv["r"],
      _["P"] = P);
  }
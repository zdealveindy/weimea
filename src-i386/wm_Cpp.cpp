#include <RcppArmadillo.h>
#include <Rmath.h>
// [[Rcpp::depends("RcppArmadillo")]]

using namespace Rcpp;
using namespace arma;



// // [[Rcpp::export()]]
// NumericMatrix submat_rcpp(NumericMatrix X, LogicalVector condition) {
//   int n=X.nrow(), k=X.ncol();
//   NumericMatrix out(n,sum (condition));
//   for (int i = 0, j = 0; i < k; i++) {
//     if(condition[i]) {
//       out(_,j) = X(_,i);
//       j = j+1;
//     }
//   }
//   return(out);
// }

// //[[Rcpp::export()]]
// NumericVector subvec_rcpp (NumericVector x, LogicalVector condition) {
//   int n = x.size ();
//   NumericVector out (sum (condition));
//   for (int i = 0, j = 0; i < n; i++) {
//     if (condition [i]) {
//       out(j) = x (i);
//       j = j+1;
//     }
//   }
//   return (out);
// }


// //[[Rcpp::export()]]
// NumericVector rowsum_rcpp (NumericMatrix X) {
//   NumericVector rowsum (X.nrow ());
//   for (int i = 0; i < X.nrow (); i++){
//     double total = 0;
//     for (int j = 0; j < X.ncol (); j++){
//       total += X (i,j);
//     }
//     rowsum[i] = total;
//   }
//   return (rowsum);
// }

// // [[Rcpp::export()]]
// NumericMatrix stand_tot_rcpp (NumericMatrix sitspe) {
//   
//   int sitspe_nrow = sitspe.nrow (), sitspe_ncol = sitspe.ncol ();
//   
//   NumericVector rowsum_sitspe = rowsum_rcpp (sitspe);
//   for (int i = 0; i < sitspe_nrow; i++){
//     for (int j = 0; j < sitspe_ncol; j++){
//       sitspe(i,j) = sitspe(i,j)/rowsum_sitspe[i];
//     }
//   }
//   return (sitspe);
// }

// //[[Rcpp::export()]]
// NumericMatrix wm_Cpp (NumericMatrix sitspe, NumericMatrix speatt) {
//   int sitspe_nrow = sitspe.nrow (), speatt_ncol = speatt.ncol ();
//   
//   NumericMatrix out (sitspe_nrow, speatt_ncol);
//   for (int sa = 0; sa < speatt_ncol; sa++) {
//     NumericMatrix sitspe_temp = stand_tot_rcpp (submat_rcpp (sitspe, !is_na (speatt(_, sa))));
//     NumericVector speatt_temp = subvec_rcpp (speatt (_, sa), !is_na (speatt (_, sa)));
//     NumericMatrix sitspe_speatt (sitspe.nrow (), sitspe.ncol ());
//     for (int ro = 0; ro < sitspe_temp.nrow (); ro++){
//       for (int co = 0; co < sitspe_temp.ncol (); co++) {
//         sitspe_speatt (ro, co) = sitspe_temp (ro,co)*speatt_temp (co);
//       }
//     }
//     out (_, sa) = rowsum_rcpp (sitspe_speatt);
//   }
//   return (out);
// }

//[[Rcpp::export()]]
arma::mat stand_tot (arma::mat sitspe) {
  for (int i = 0; i < sitspe.n_rows; i++){
    sitspe.row(i) /= sum (sitspe.row (i));
  }
  return (sitspe);
}

//[[Rcpp::export()]]
arma::mat wm_rcpp (arma::mat sitspe, arma::mat speatt) {
  int sitspe_nrow = sitspe.n_rows, speatt_ncol = speatt.n_cols;
  arma::mat sitspe_speatt (sitspe_nrow, speatt_ncol);
  for (int sa = 0; sa < speatt_ncol; sa++) {
    arma::mat sitspe_temp = sitspe.cols (find_finite (speatt.col(sa)));
    arma::vec speatt_temp = speatt.col (sa);
    speatt_temp = speatt_temp.elem (find_finite (speatt_temp));
    sitspe_speatt.col (sa) =  stand_tot (sitspe_temp) * speatt_temp;
  }
  return sitspe_speatt;
}

// [[Rcpp::export()]]
int count_if(LogicalVector x) {
  int counter = 0;
  for(int i = 0; i < x.size(); i++) {
    if(x[i] == TRUE) {
      counter++;
    }
  }
  return counter;
}

//[[Rcpp::export()]]
List test_LR_cor (arma::mat sitspe, arma::mat speatt, arma::mat env, CharacterVector cor_coef, double perm) {
  arma::mat sitspe_temp = sitspe (find_finite (env), find_finite (speatt));
  arma::mat env_temp = env.elem (find_finite (env));
  arma::vec M_art = wm_rcpp (sitspe_temp, wm_rcpp (sitspe_temp.t (), env_temp));
  double r_obs = as_scalar (cor (M_art, env_temp));  // calculates Pearson's r correlation coefficient
  double t_obs = r_obs*sqrt ((env_temp.size () - 2)/(1 - pow (r_obs, 2.0)));  // calculates t-value according to Student's t formula
  arma::vec r_exp (perm+1);
  arma::vec t_exp (perm+1);
  for (int nperm = 0; nperm < perm; nperm++){
    arma::mat env_rand = shuffle (env_temp);
    arma::vec M_exp = wm_rcpp (sitspe_temp, wm_rcpp (sitspe_temp.t(), env_rand));
    r_exp [nperm] = as_scalar (cor (M_exp, env_rand));
    t_exp [nperm] = r_exp [nperm]*sqrt ((env_rand.size () - 2)/(1 - pow (r_exp[nperm], 2)));
  }
  r_exp (perm) = r_obs; //observed values is put on the last place of the vector
  t_exp (perm) = t_obs; //observed values is put on the last place of the vector
  double P = sum (abs (t_exp) >= abs (t_obs))/(perm+1);
  return List::create (
      _["type"] = "cor",
      _["no.samples"] = sitspe.n_rows,
      _["no.species"] = sitspe.n_cols,
      _["cor.coef"] = cor_coef,
      _["statistic"] = t_obs,
      _["estimate"] = r_obs,
      _["cor"] = P
  );
}


//[[Rcpp::export()]]
bool is_in (CharacterVector x, CharacterVector table){
  bool is_in_temp = FALSE;
  for (int i = 0; i < x.size (); i++){
    for (int j = 0; j < table.size (); j++){
      if (x(i) == table(j)) is_in_temp = TRUE;
      };
  };
  return is_in_temp;
}

  
//[[Rcpp::export()]]
List test_MR_cor (arma::mat sitspe, arma::mat speatt, arma::mat env, CharacterVector test, CharacterVector cor_coef, double perm, double testLR_P, double testLR_perm) {
  arma::mat sitspe_temp = sitspe.submat (find_finite (env), find_finite (speatt));
  arma::mat speatt_temp = speatt.elem (find_finite (speatt));
  arma::mat env_temp = env.elem (find_finite (env));
  arma::vec M = wm_rcpp (sitspe_temp, speatt_temp);
  double no_samples = sitspe_temp.n_rows;
  double r_obs = as_scalar (cor (M, env_temp));  // calculates Pearson's r correlation coefficient
  double t_obs = r_obs*sqrt ((env_temp.size () - 2)/(1 - pow (r_obs, 2.0)));  // calculates t-value according to Student's t formula
  double P_par = R::pt (t_obs, no_samples-2, TRUE, FALSE);
  P_par = 2*(min (NumericVector::create (P_par, 1-P_par)));
  double P_sta = NA_REAL;
  double P_mod = NA_REAL;
  double P_LR = NA_REAL;
  double P_two = NA_REAL;
  
  if (is_in (CharacterVector::create ("standard", "twostep"), test)){
    arma::vec r_exp_sta (perm+1);
    arma::vec t_exp_sta (perm+1);
    for (int nperm = 0; nperm < perm; nperm++){
      arma::mat env_rand = shuffle (env_temp);
      r_exp_sta [nperm] = as_scalar (cor (M, env_rand));  // original M with randomized R
      t_exp_sta [nperm] = r_exp_sta [nperm]*sqrt ((env_rand.size () - 2)/(1 - pow (r_exp_sta[nperm], 2)));
    }
    r_exp_sta (perm) = r_obs; //observed values is put on the last place of the vector
    t_exp_sta (perm) = t_obs; //observed values is put on the last place of the vector
    P_sta = sum (abs (t_exp_sta) >= abs (t_obs))/(perm+1);
  };
  
  if (is_in (CharacterVector::create ("modified", "twostep"), test)){
    arma::vec r_exp_mod (perm+1);
    arma::vec t_exp_mod (perm+1);
    for (int nperm = 0; nperm < perm; nperm++){
      arma::vec M_rand = wm_rcpp (sitspe_temp, shuffle (speatt_temp));
      r_exp_mod [nperm] = as_scalar (cor (M_rand, env_temp));  // original M with randomized R
      t_exp_mod [nperm] = r_exp_mod [nperm]*sqrt ((env_temp.size () - 2)/(1 - pow (r_exp_mod[nperm], 2)));
    }
    r_exp_mod (perm) = r_obs; //observed values is put on the last place of the vector
    t_exp_mod (perm) = t_obs; //observed values is put on the last place of the vector
    P_mod = sum (abs (t_exp_mod) >= abs (t_obs))/(perm+1);
  };

  if (is_in ("twostep", test)){
    arma::vec M_art = wm_rcpp (sitspe_temp, wm_rcpp (sitspe_temp.t (), env_temp));
    double r_obs_LR = as_scalar (cor (M_art, env_temp));  // calculates Pearson's r correlation coefficient
    double t_obs_LR = r_obs_LR*sqrt ((env_temp.size () - 2)/(1 - pow (r_obs_LR, 2.0)));  // calculates t-value according to Student's t formula
    arma::vec r_exp_LR (testLR_perm+1);
    arma::vec t_exp_LR (testLR_perm+1);
    for (int nperm = 0; nperm < testLR_perm; nperm++){
      arma::mat env_rand = shuffle (env_temp);
      arma::vec M_exp = wm_rcpp (sitspe_temp, wm_rcpp (sitspe_temp.t(), env_rand));
      r_exp_LR [nperm] = as_scalar (cor (M_exp, env_rand));
      t_exp_LR [nperm] = r_exp_LR [nperm]*sqrt ((env_rand.size () - 2)/(1 - pow (r_exp_LR[nperm], 2)));
    }
    r_exp_LR (perm) = r_obs_LR; //observed values is put on the last place of the vector
    t_exp_LR (perm) = t_obs_LR; //observed values is put on the last place of the vector
    P_LR = sum (abs (t_exp_LR) >= abs (t_obs_LR))/(perm+1);
    if (P_LR <= testLR_P) {P_two = P_mod;} else {P_two = P_sta;};
  };
  
  return List::create (
    _["r.obs"] = r_obs,
    _["t.obs"] = t_obs,
    _["P.par"] = P_par,
    _["P.sta"] = P_sta,
    _["P.mod"] = P_mod,
    _["P.LR"] = P_LR,
    _["P.two"] = P_two
  );
}


// Extended version of fastLm from RcppArmadillo
// [[Rcpp::export]]
List fastLm_wm (const arma::mat& X, const arma::colvec& y, bool intercept_included) {
  int n = X.n_rows, k = X.n_cols;
  
  arma::colvec coef = arma::solve(X, y);    // fit model y ~ X
  arma::colvec res  = y - X*coef;           // residuals
  arma::colvec fit = X*coef;                // fitted values
  double mss = sum (pow (fit - mean (fit), 2.0));    // sum of squares   
  double rss = sum (pow (res - 0, 2.0));              // residual sum of squares  NEED TO SOLVE THE PROBLEMS WITH DOUBLE vs MAT conversion here!
  double rsq = mss/(mss + rss);           // r2
  int mdf;
  if (intercept_included) {mdf = k-1;} else {mdf = k;};
  int rdf = n - k;      // residual number of degrees of freedom, should be n - no.vars - 1 (but X contains 1 column with 1's for intercept
  double F = (mss/mdf)/(rss/rdf);
//   // std.errors of coefficients
//   double s2 = std::inner_product(res.begin(), res.end(), res.begin(), 0.0)/(n - k);
//   
//   arma::colvec std_err = arma::sqrt(s2 * arma::diagvec(arma::pinv(arma::trans(X)*X)));  
  return List::create(_["coefficients"] = coef,
  // not needed       _["stderr"]       = std_err,
                      _["df.residual"]  = rdf,
                      _["df.model"] = mdf,
  // not needed       _["residuals"] = res,
                      _["F.value"] = F,
                      _["R.squared"] = rsq
                      );
}


//[[Rcpp::export()]]
List test_LR_lm (arma::mat sitspe, arma::mat speatt, arma::mat env, double perm) {
  arma::mat sitspe_temp = sitspe (find_finite (env), find_finite (speatt));
  arma::mat env_temp = env.elem (find_finite (env));
  arma::vec M_art = wm_rcpp (sitspe_temp, wm_rcpp (sitspe_temp.t (), env_temp));
  List lm_obs = fastLm_wm (join_rows (arma::ones (env_temp.n_rows), env_temp), M_art, TRUE);
  double rsq_obs = lm_obs["R.squared"];
  double F_obs = lm_obs["F.value"];

  arma::vec rsq_exp(perm+1);
  arma::vec F_exp (perm+1);
  for (int nperm = 0; nperm < perm; nperm++){
    arma::mat env_rand = shuffle (env_temp);
    arma::vec M_exp = wm_rcpp (sitspe_temp, wm_rcpp (sitspe_temp.t(), env_rand));
    List lm_exp = fastLm_wm (join_rows (arma::ones (env_temp.n_rows), env_rand), M_exp, TRUE);
    rsq_exp (nperm) = lm_exp["R.squared"];
    F_exp (nperm) = lm_exp["F.value"];
  }
  rsq_exp (perm) = rsq_obs; //observed values is put on the last place of the vector
  F_exp (perm) = F_obs; //observed values is put on the last place of the vector
  double P = sum (abs (F_exp) >= abs (F_obs))/(perm+1);
  return List::create (
      _["type"] = "lm",
      _["no.samples"] = sitspe.n_rows,
      _["no.species"] = sitspe.n_cols,
      _["reg.coef"] = lm_obs["coefficients"],
      _["statistic"] = F_obs,
      _["estimate"] = rsq_obs,
      _["lm"] = P
  );
}

// [[Rcpp::export]]

arma::uvec keep_rows (arma::mat X){
int n = X.n_rows,k = X.n_cols;
// create keep vector
arma::vec keep = arma::ones<arma::vec>(n);
for (int j = 0; j < k; j++) 
  for (int i = 0; i < n; i++) 
    if (keep[i] && !arma::is_finite(X(i,j))) keep[i] = 0;
    return (find (keep == 1));
}

//[[Rcpp::export()]]
List test_MR_lm (arma::mat sitspe, arma::mat speatt, arma::mat env, CharacterVector test, const char* dependence, double perm, double testLR_P, double testLR_perm) {
  arma::mat sitspe_temp = sitspe.submat (keep_rows (env), keep_rows (speatt));
  arma::mat speatt_temp = speatt.rows (keep_rows (speatt));
  arma::mat env_temp = env.rows (keep_rows (env));
  arma::mat M = wm_rcpp (sitspe_temp, speatt_temp);
  List lm_obs;
 if (!strncmp (dependence, "M ~ env", 1)) {
   lm_obs = fastLm_wm (join_rows (arma::ones (env_temp.n_rows), env_temp), M, TRUE);
   }
   else {
    lm_obs = fastLm_wm (join_rows (arma::ones (M.n_rows), M), env_temp, TRUE);    
  }

  double rsq_obs = lm_obs["R.squared"];
  double rsq_adj_sta = NA_REAL;
  double rsq_adj_mod = NA_REAL;
  double F_obs = lm_obs["F.value"];
  double P_par = R::pf (F_obs, lm_obs["df.model"], lm_obs["df.residual"], FALSE, FALSE);
  double P_sta = NA_REAL;
  double P_mod = NA_REAL;
  double P_LR = NA_REAL;
  double P_two = NA_REAL;
 
  if (is_in (CharacterVector::create ("standard", "twostep"), test)){
    arma::vec rsq_exp_sta (perm);
    arma::vec F_exp_sta (perm+1);
    for (int nperm = 0; nperm < perm; nperm++){
      arma::mat env_rand = shuffle (env_temp, 0);  // rows of env_temp are shuffled
      List lm_exp;
      if (!strncmp (dependence, "M ~ env", 1)) {
        lm_exp = fastLm_wm (join_rows (arma::ones (env_rand.n_rows), env_rand), M, TRUE);
      } 
      else {
        lm_exp = fastLm_wm (join_rows (arma::ones (M.n_rows), M), env_rand, TRUE);
      }
      rsq_exp_sta (nperm) = lm_exp["R.squared"];  // r.sq for calculation of permutation-based adjusted r.sq
      F_exp_sta (nperm) = lm_exp["F.value"];
    }
    F_exp_sta (perm) = F_obs; //observed values is put on the last place of the vector
    P_sta = sum (abs (F_exp_sta) >= abs (F_obs))/(perm+1);
    rsq_adj_sta = 1 - (1/(1-mean (rsq_exp_sta))*(1-rsq_obs));
  };
  
  if (is_in (CharacterVector::create ("modified", "twostep"), test)){
    arma::vec rsq_exp_mod (perm);  // for adjusted R2
    arma::vec F_exp_mod (perm+1);
    for (int nperm = 0; nperm < perm; nperm++){
      arma::mat M_rand = wm_rcpp (sitspe_temp, shuffle (speatt_temp, 0));
      List lm_exp;
      if (!strncmp (dependence, "M ~ env", 1)) {  
        lm_exp = fastLm_wm (join_rows (arma::ones (env_temp.n_rows), env_temp), M_rand, TRUE);
      }
      else {
        lm_exp = fastLm_wm (join_rows (arma::ones (M_rand.n_rows), M_rand), env_temp, TRUE);
      }
      rsq_exp_mod [nperm] = lm_exp["R.squared"];  // r.sq for calculation of permutation-based r.sq
      F_exp_mod [nperm] = lm_exp["F.value"];
    }
    F_exp_mod (perm) = F_obs; //observed values is put on the last place of the vector
    P_mod = sum (abs (F_exp_mod) >= abs (F_obs))/(perm+1);
    rsq_adj_mod = 1 - (1/(1-mean (rsq_exp_mod))*(1-rsq_obs));
  };
  
  if (is_in ("twostep", test)){
    arma::mat M_art = wm_rcpp (sitspe_temp, wm_rcpp (sitspe_temp.t (), env_temp));
    List lm_obs_LR;
    if (!strncmp (dependence, "M ~ env", 1)) {
      lm_obs_LR = fastLm_wm (join_rows (arma::ones (env_temp.n_rows), env_temp), M_art, TRUE);
    }
    else {
      lm_obs_LR = fastLm_wm (join_rows (arma::ones (M_art.n_rows), M_art), env_temp, TRUE);
    }
    double rsq_obs_LR = lm_obs_LR["R.squared"];
    double F_obs_LR = lm_obs_LR["F.value"];
    arma::vec rsq_exp_LR (testLR_perm+1);
    arma::vec F_exp_LR (testLR_perm+1);
    for (int nperm = 0; nperm < testLR_perm; nperm++){
      arma::mat env_rand = shuffle (env_temp);
      arma::vec M_exp = wm_rcpp (sitspe_temp, wm_rcpp (sitspe_temp.t(), env_rand));
      List lm_exp_LR = fastLm_wm (join_rows (arma::ones (env_rand.n_rows), env_rand), M_exp, TRUE);
      rsq_exp_LR [nperm] = lm_exp_LR["R.squared"];
      F_exp_LR [nperm] = lm_exp_LR["F.value"];
    }
    rsq_exp_LR (perm) = rsq_obs_LR; //observed values is put on the last place of the vector
    F_exp_LR (perm) = F_obs_LR; //observed values is put on the last place of the vector
    P_LR = sum (abs (F_exp_LR) >= abs (F_obs_LR))/(perm+1);
    if (P_LR <= testLR_P) {P_two = P_mod;} else {P_two = P_sta;};
  };
  
  return List::create (
      _["coef"] = lm_obs["coefficients"],
      _["rsq.obs"] = rsq_obs,
      _["rsq.adj.sta"] = rsq_adj_sta,
      _["rsq.adj.mod"] = rsq_adj_mod,
      _["F.obs"] = F_obs,
      _["P.par"] = P_par,
      _["P.sta"] = P_sta,
      _["P.mod"] = P_mod,
      _["P.LR"] = P_LR,
      _["P.two"] = P_two
);
}

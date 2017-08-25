// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// stand_tot
arma::mat stand_tot(arma::mat sitspe);
RcppExport SEXP _weimea_stand_tot(SEXP sitspeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type sitspe(sitspeSEXP);
    rcpp_result_gen = Rcpp::wrap(stand_tot(sitspe));
    return rcpp_result_gen;
END_RCPP
}
// wm_rcpp
arma::mat wm_rcpp(arma::mat sitspe, arma::mat speatt);
RcppExport SEXP _weimea_wm_rcpp(SEXP sitspeSEXP, SEXP speattSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type sitspe(sitspeSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type speatt(speattSEXP);
    rcpp_result_gen = Rcpp::wrap(wm_rcpp(sitspe, speatt));
    return rcpp_result_gen;
END_RCPP
}
// count_if
int count_if(LogicalVector x);
RcppExport SEXP _weimea_count_if(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< LogicalVector >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(count_if(x));
    return rcpp_result_gen;
END_RCPP
}
// test_LR_cor
List test_LR_cor(arma::mat sitspe, arma::mat speatt, arma::mat env, CharacterVector cor_coef, double perm);
RcppExport SEXP _weimea_test_LR_cor(SEXP sitspeSEXP, SEXP speattSEXP, SEXP envSEXP, SEXP cor_coefSEXP, SEXP permSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type sitspe(sitspeSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type speatt(speattSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type env(envSEXP);
    Rcpp::traits::input_parameter< CharacterVector >::type cor_coef(cor_coefSEXP);
    Rcpp::traits::input_parameter< double >::type perm(permSEXP);
    rcpp_result_gen = Rcpp::wrap(test_LR_cor(sitspe, speatt, env, cor_coef, perm));
    return rcpp_result_gen;
END_RCPP
}
// is_in
bool is_in(CharacterVector x, CharacterVector table);
RcppExport SEXP _weimea_is_in(SEXP xSEXP, SEXP tableSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< CharacterVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< CharacterVector >::type table(tableSEXP);
    rcpp_result_gen = Rcpp::wrap(is_in(x, table));
    return rcpp_result_gen;
END_RCPP
}
// test_MR_cor_pear
List test_MR_cor_pear(arma::mat sitspe, arma::mat speatt, arma::mat env, CharacterVector test, double perm, double testLR_P, double testLR_perm);
RcppExport SEXP _weimea_test_MR_cor_pear(SEXP sitspeSEXP, SEXP speattSEXP, SEXP envSEXP, SEXP testSEXP, SEXP permSEXP, SEXP testLR_PSEXP, SEXP testLR_permSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type sitspe(sitspeSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type speatt(speattSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type env(envSEXP);
    Rcpp::traits::input_parameter< CharacterVector >::type test(testSEXP);
    Rcpp::traits::input_parameter< double >::type perm(permSEXP);
    Rcpp::traits::input_parameter< double >::type testLR_P(testLR_PSEXP);
    Rcpp::traits::input_parameter< double >::type testLR_perm(testLR_permSEXP);
    rcpp_result_gen = Rcpp::wrap(test_MR_cor_pear(sitspe, speatt, env, test, perm, testLR_P, testLR_perm));
    return rcpp_result_gen;
END_RCPP
}
// fastLm_wm
List fastLm_wm(const arma::mat& X, const arma::colvec& y, bool intercept_included);
RcppExport SEXP _weimea_fastLm_wm(SEXP XSEXP, SEXP ySEXP, SEXP intercept_includedSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< bool >::type intercept_included(intercept_includedSEXP);
    rcpp_result_gen = Rcpp::wrap(fastLm_wm(X, y, intercept_included));
    return rcpp_result_gen;
END_RCPP
}
// test_LR_lm
List test_LR_lm(arma::mat sitspe, arma::mat speatt, arma::mat env, double perm);
RcppExport SEXP _weimea_test_LR_lm(SEXP sitspeSEXP, SEXP speattSEXP, SEXP envSEXP, SEXP permSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type sitspe(sitspeSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type speatt(speattSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type env(envSEXP);
    Rcpp::traits::input_parameter< double >::type perm(permSEXP);
    rcpp_result_gen = Rcpp::wrap(test_LR_lm(sitspe, speatt, env, perm));
    return rcpp_result_gen;
END_RCPP
}
// keep_rows
arma::uvec keep_rows(arma::mat X);
RcppExport SEXP _weimea_keep_rows(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(keep_rows(X));
    return rcpp_result_gen;
END_RCPP
}
// test_MR_lm
List test_MR_lm(arma::mat sitspe, arma::mat speatt, arma::mat env, CharacterVector test, const char* dependence, double perm, double testLR_P, double testLR_perm);
RcppExport SEXP _weimea_test_MR_lm(SEXP sitspeSEXP, SEXP speattSEXP, SEXP envSEXP, SEXP testSEXP, SEXP dependenceSEXP, SEXP permSEXP, SEXP testLR_PSEXP, SEXP testLR_permSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type sitspe(sitspeSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type speatt(speattSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type env(envSEXP);
    Rcpp::traits::input_parameter< CharacterVector >::type test(testSEXP);
    Rcpp::traits::input_parameter< const char* >::type dependence(dependenceSEXP);
    Rcpp::traits::input_parameter< double >::type perm(permSEXP);
    Rcpp::traits::input_parameter< double >::type testLR_P(testLR_PSEXP);
    Rcpp::traits::input_parameter< double >::type testLR_perm(testLR_permSEXP);
    rcpp_result_gen = Rcpp::wrap(test_MR_lm(sitspe, speatt, env, test, dependence, perm, testLR_P, testLR_perm));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_weimea_stand_tot", (DL_FUNC) &_weimea_stand_tot, 1},
    {"_weimea_wm_rcpp", (DL_FUNC) &_weimea_wm_rcpp, 2},
    {"_weimea_count_if", (DL_FUNC) &_weimea_count_if, 1},
    {"_weimea_test_LR_cor", (DL_FUNC) &_weimea_test_LR_cor, 5},
    {"_weimea_is_in", (DL_FUNC) &_weimea_is_in, 2},
    {"_weimea_test_MR_cor_pear", (DL_FUNC) &_weimea_test_MR_cor_pear, 7},
    {"_weimea_fastLm_wm", (DL_FUNC) &_weimea_fastLm_wm, 3},
    {"_weimea_test_LR_lm", (DL_FUNC) &_weimea_test_LR_lm, 4},
    {"_weimea_keep_rows", (DL_FUNC) &_weimea_keep_rows, 1},
    {"_weimea_test_MR_lm", (DL_FUNC) &_weimea_test_MR_lm, 8},
    {NULL, NULL, 0}
};

RcppExport void R_init_weimea(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
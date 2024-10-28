// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// rnkCalc
NumericMatrix rnkCalc(NumericMatrix logrhonk, NumericVector lse, double N, double K);
RcppExport SEXP _VICatMix_rnkCalc(SEXP logrhonkSEXP, SEXP lseSEXP, SEXP NSEXP, SEXP KSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type logrhonk(logrhonkSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type lse(lseSEXP);
    Rcpp::traits::input_parameter< double >::type N(NSEXP);
    Rcpp::traits::input_parameter< double >::type K(KSEXP);
    rcpp_result_gen = Rcpp::wrap(rnkCalc(logrhonk, lse, N, K));
    return rcpp_result_gen;
END_RCPP
}
// ElogphiCalc
arma::cube ElogphiCalc(arma::cube eps, double K, double D, double N, double maxNCat, arma::mat X);
RcppExport SEXP _VICatMix_ElogphiCalc(SEXP epsSEXP, SEXP KSEXP, SEXP DSEXP, SEXP NSEXP, SEXP maxNCatSEXP, SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::cube >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< double >::type K(KSEXP);
    Rcpp::traits::input_parameter< double >::type D(DSEXP);
    Rcpp::traits::input_parameter< double >::type N(NSEXP);
    Rcpp::traits::input_parameter< double >::type maxNCat(maxNCatSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(ElogphiCalc(eps, K, D, N, maxNCat, X));
    return rcpp_result_gen;
END_RCPP
}
// ElogphiLCalc
arma::cube ElogphiLCalc(arma::cube eps, double K, double D, double maxNCat);
RcppExport SEXP _VICatMix_ElogphiLCalc(SEXP epsSEXP, SEXP KSEXP, SEXP DSEXP, SEXP maxNCatSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::cube >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< double >::type K(KSEXP);
    Rcpp::traits::input_parameter< double >::type D(DSEXP);
    Rcpp::traits::input_parameter< double >::type maxNCat(maxNCatSEXP);
    rcpp_result_gen = Rcpp::wrap(ElogphiLCalc(eps, K, D, maxNCat));
    return rcpp_result_gen;
END_RCPP
}
// logrhonkCalc
arma::mat logrhonkCalc(arma::vec Elogpi, arma::cube Elogphi, double K, double D, double N);
RcppExport SEXP _VICatMix_logrhonkCalc(SEXP ElogpiSEXP, SEXP ElogphiSEXP, SEXP KSEXP, SEXP DSEXP, SEXP NSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type Elogpi(ElogpiSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type Elogphi(ElogphiSEXP);
    Rcpp::traits::input_parameter< double >::type K(KSEXP);
    Rcpp::traits::input_parameter< double >::type D(DSEXP);
    Rcpp::traits::input_parameter< double >::type N(NSEXP);
    rcpp_result_gen = Rcpp::wrap(logrhonkCalc(Elogpi, Elogphi, K, D, N));
    return rcpp_result_gen;
END_RCPP
}
// epsCalc
arma::cube epsCalc(double K, double maxNCat, double D, double N, arma::mat prioreps, arma::mat X, arma::mat rnk);
RcppExport SEXP _VICatMix_epsCalc(SEXP KSEXP, SEXP maxNCatSEXP, SEXP DSEXP, SEXP NSEXP, SEXP priorepsSEXP, SEXP XSEXP, SEXP rnkSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type K(KSEXP);
    Rcpp::traits::input_parameter< double >::type maxNCat(maxNCatSEXP);
    Rcpp::traits::input_parameter< double >::type D(DSEXP);
    Rcpp::traits::input_parameter< double >::type N(NSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type prioreps(priorepsSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type rnk(rnkSEXP);
    rcpp_result_gen = Rcpp::wrap(epsCalc(K, maxNCat, D, N, prioreps, X, rnk));
    return rcpp_result_gen;
END_RCPP
}
// firstepsCalc
arma::cube firstepsCalc(double K, double maxNCat, double D, double N, arma::mat prioreps, arma::mat X, arma::vec clusterInit);
RcppExport SEXP _VICatMix_firstepsCalc(SEXP KSEXP, SEXP maxNCatSEXP, SEXP DSEXP, SEXP NSEXP, SEXP priorepsSEXP, SEXP XSEXP, SEXP clusterInitSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type K(KSEXP);
    Rcpp::traits::input_parameter< double >::type maxNCat(maxNCatSEXP);
    Rcpp::traits::input_parameter< double >::type D(DSEXP);
    Rcpp::traits::input_parameter< double >::type N(NSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type prioreps(priorepsSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type clusterInit(clusterInitSEXP);
    rcpp_result_gen = Rcpp::wrap(firstepsCalc(K, maxNCat, D, N, prioreps, X, clusterInit));
    return rcpp_result_gen;
END_RCPP
}
// CpriorepsCalc
arma::mat CpriorepsCalc(arma::mat prioreps, double K, double D, arma::vec nCat);
RcppExport SEXP _VICatMix_CpriorepsCalc(SEXP priorepsSEXP, SEXP KSEXP, SEXP DSEXP, SEXP nCatSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type prioreps(priorepsSEXP);
    Rcpp::traits::input_parameter< double >::type K(KSEXP);
    Rcpp::traits::input_parameter< double >::type D(DSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type nCat(nCatSEXP);
    rcpp_result_gen = Rcpp::wrap(CpriorepsCalc(prioreps, K, D, nCat));
    return rcpp_result_gen;
END_RCPP
}
// CpostepsCalc
arma::mat CpostepsCalc(arma::cube eps, double K, double D, double maxNCat);
RcppExport SEXP _VICatMix_CpostepsCalc(SEXP epsSEXP, SEXP KSEXP, SEXP DSEXP, SEXP maxNCatSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::cube >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< double >::type K(KSEXP);
    Rcpp::traits::input_parameter< double >::type D(DSEXP);
    Rcpp::traits::input_parameter< double >::type maxNCat(maxNCatSEXP);
    rcpp_result_gen = Rcpp::wrap(CpostepsCalc(eps, K, D, maxNCat));
    return rcpp_result_gen;
END_RCPP
}
// sumDElogphiCalc
arma::mat sumDElogphiCalc(arma::cube Elogphi, double K, double D, double N);
RcppExport SEXP _VICatMix_sumDElogphiCalc(SEXP ElogphiSEXP, SEXP KSEXP, SEXP DSEXP, SEXP NSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::cube >::type Elogphi(ElogphiSEXP);
    Rcpp::traits::input_parameter< double >::type K(KSEXP);
    Rcpp::traits::input_parameter< double >::type D(DSEXP);
    Rcpp::traits::input_parameter< double >::type N(NSEXP);
    rcpp_result_gen = Rcpp::wrap(sumDElogphiCalc(Elogphi, K, D, N));
    return rcpp_result_gen;
END_RCPP
}
// priorepsminusoneCalc
arma::cube priorepsminusoneCalc(arma::mat prioreps, double K, double D, double maxNCat);
RcppExport SEXP _VICatMix_priorepsminusoneCalc(SEXP priorepsSEXP, SEXP KSEXP, SEXP DSEXP, SEXP maxNCatSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type prioreps(priorepsSEXP);
    Rcpp::traits::input_parameter< double >::type K(KSEXP);
    Rcpp::traits::input_parameter< double >::type D(DSEXP);
    Rcpp::traits::input_parameter< double >::type maxNCat(maxNCatSEXP);
    rcpp_result_gen = Rcpp::wrap(priorepsminusoneCalc(prioreps, K, D, maxNCat));
    return rcpp_result_gen;
END_RCPP
}
// epsminusoneCalc
arma::cube epsminusoneCalc(arma::cube eps, double K, double D, double maxNCat);
RcppExport SEXP _VICatMix_epsminusoneCalc(SEXP epsSEXP, SEXP KSEXP, SEXP DSEXP, SEXP maxNCatSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::cube >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< double >::type K(KSEXP);
    Rcpp::traits::input_parameter< double >::type D(DSEXP);
    Rcpp::traits::input_parameter< double >::type maxNCat(maxNCatSEXP);
    rcpp_result_gen = Rcpp::wrap(epsminusoneCalc(eps, K, D, maxNCat));
    return rcpp_result_gen;
END_RCPP
}
// epsminuspriorepsCalc
arma::cube epsminuspriorepsCalc(arma::cube eps, arma::mat prioreps, double K, double D, double maxNCat);
RcppExport SEXP _VICatMix_epsminuspriorepsCalc(SEXP epsSEXP, SEXP priorepsSEXP, SEXP KSEXP, SEXP DSEXP, SEXP maxNCatSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::cube >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type prioreps(priorepsSEXP);
    Rcpp::traits::input_parameter< double >::type K(KSEXP);
    Rcpp::traits::input_parameter< double >::type D(DSEXP);
    Rcpp::traits::input_parameter< double >::type maxNCat(maxNCatSEXP);
    rcpp_result_gen = Rcpp::wrap(epsminuspriorepsCalc(eps, prioreps, K, D, maxNCat));
    return rcpp_result_gen;
END_RCPP
}
// cmatrixCalc
NumericMatrix cmatrixCalc(NumericMatrix nullphi, NumericMatrix X, NumericVector c, double N, double D);
RcppExport SEXP _VICatMix_cmatrixCalc(SEXP nullphiSEXP, SEXP XSEXP, SEXP cSEXP, SEXP NSEXP, SEXP DSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type nullphi(nullphiSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type X(XSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type c(cSEXP);
    Rcpp::traits::input_parameter< double >::type N(NSEXP);
    Rcpp::traits::input_parameter< double >::type D(DSEXP);
    rcpp_result_gen = Rcpp::wrap(cmatrixCalc(nullphi, X, c, N, D));
    return rcpp_result_gen;
END_RCPP
}
// ElogthetaCalcCat
NumericMatrix ElogthetaCalcCat(NumericMatrix beta, double K, double J);
RcppExport SEXP _VICatMix_ElogthetaCalcCat(SEXP betaSEXP, SEXP KSEXP, SEXP JSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< double >::type K(KSEXP);
    Rcpp::traits::input_parameter< double >::type J(JSEXP);
    rcpp_result_gen = Rcpp::wrap(ElogthetaCalcCat(beta, K, J));
    return rcpp_result_gen;
END_RCPP
}
// betaCalc
NumericMatrix betaCalc(NumericVector priorbeta, NumericVector y, double K, double J, double N, NumericMatrix rnk);
RcppExport SEXP _VICatMix_betaCalc(SEXP priorbetaSEXP, SEXP ySEXP, SEXP KSEXP, SEXP JSEXP, SEXP NSEXP, SEXP rnkSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type priorbeta(priorbetaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< double >::type K(KSEXP);
    Rcpp::traits::input_parameter< double >::type J(JSEXP);
    Rcpp::traits::input_parameter< double >::type N(NSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type rnk(rnkSEXP);
    rcpp_result_gen = Rcpp::wrap(betaCalc(priorbeta, y, K, J, N, rnk));
    return rcpp_result_gen;
END_RCPP
}
// CpriorbetaCalc
NumericVector CpriorbetaCalc(NumericVector priorbeta, double K, double J);
RcppExport SEXP _VICatMix_CpriorbetaCalc(SEXP priorbetaSEXP, SEXP KSEXP, SEXP JSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type priorbeta(priorbetaSEXP);
    Rcpp::traits::input_parameter< double >::type K(KSEXP);
    Rcpp::traits::input_parameter< double >::type J(JSEXP);
    rcpp_result_gen = Rcpp::wrap(CpriorbetaCalc(priorbeta, K, J));
    return rcpp_result_gen;
END_RCPP
}
// CpostbetaCalc
NumericVector CpostbetaCalc(NumericMatrix beta, double K, double J);
RcppExport SEXP _VICatMix_CpostbetaCalc(SEXP betaSEXP, SEXP KSEXP, SEXP JSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< double >::type K(KSEXP);
    Rcpp::traits::input_parameter< double >::type J(JSEXP);
    rcpp_result_gen = Rcpp::wrap(CpostbetaCalc(beta, K, J));
    return rcpp_result_gen;
END_RCPP
}
// respthetaCalc
NumericMatrix respthetaCalc(NumericMatrix Elogtheta, NumericMatrix rnk, NumericVector y, double N, double K);
RcppExport SEXP _VICatMix_respthetaCalc(SEXP ElogthetaSEXP, SEXP rnkSEXP, SEXP ySEXP, SEXP NSEXP, SEXP KSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type Elogtheta(ElogthetaSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type rnk(rnkSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< double >::type N(NSEXP);
    Rcpp::traits::input_parameter< double >::type K(KSEXP);
    rcpp_result_gen = Rcpp::wrap(respthetaCalc(Elogtheta, rnk, y, N, K));
    return rcpp_result_gen;
END_RCPP
}
// nullphiCalc
NumericMatrix nullphiCalc(NumericMatrix X, NumericVector nCat, double maxNCat, double D, double N);
RcppExport SEXP _VICatMix_nullphiCalc(SEXP XSEXP, SEXP nCatSEXP, SEXP maxNCatSEXP, SEXP DSEXP, SEXP NSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type X(XSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type nCat(nCatSEXP);
    Rcpp::traits::input_parameter< double >::type maxNCat(maxNCatSEXP);
    Rcpp::traits::input_parameter< double >::type D(DSEXP);
    Rcpp::traits::input_parameter< double >::type N(NSEXP);
    rcpp_result_gen = Rcpp::wrap(nullphiCalc(X, nCat, maxNCat, D, N));
    return rcpp_result_gen;
END_RCPP
}
// firstbetaCalc
NumericMatrix firstbetaCalc(NumericVector y, NumericVector priorbeta, double K, double J, double N, NumericVector clusterInit);
RcppExport SEXP _VICatMix_firstbetaCalc(SEXP ySEXP, SEXP priorbetaSEXP, SEXP KSEXP, SEXP JSEXP, SEXP NSEXP, SEXP clusterInitSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericVector >::type priorbeta(priorbetaSEXP);
    Rcpp::traits::input_parameter< double >::type K(KSEXP);
    Rcpp::traits::input_parameter< double >::type J(JSEXP);
    Rcpp::traits::input_parameter< double >::type N(NSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type clusterInit(clusterInitSEXP);
    rcpp_result_gen = Rcpp::wrap(firstbetaCalc(y, priorbeta, K, J, N, clusterInit));
    return rcpp_result_gen;
END_RCPP
}
// lognullphiCalc
arma::cube lognullphiCalc(arma::mat nullphi, arma::mat X, double K, double D, double N);
RcppExport SEXP _VICatMix_lognullphiCalc(SEXP nullphiSEXP, SEXP XSEXP, SEXP KSEXP, SEXP DSEXP, SEXP NSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type nullphi(nullphiSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< double >::type K(KSEXP);
    Rcpp::traits::input_parameter< double >::type D(DSEXP);
    Rcpp::traits::input_parameter< double >::type N(NSEXP);
    rcpp_result_gen = Rcpp::wrap(lognullphiCalc(nullphi, X, K, D, N));
    return rcpp_result_gen;
END_RCPP
}
// logrhonkCalcVarSel
arma::mat logrhonkCalcVarSel(arma::vec Elogpi, arma::cube carray, arma::mat cmatrix, double K, double D, double N);
RcppExport SEXP _VICatMix_logrhonkCalcVarSel(SEXP ElogpiSEXP, SEXP carraySEXP, SEXP cmatrixSEXP, SEXP KSEXP, SEXP DSEXP, SEXP NSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type Elogpi(ElogpiSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type carray(carraySEXP);
    Rcpp::traits::input_parameter< arma::mat >::type cmatrix(cmatrixSEXP);
    Rcpp::traits::input_parameter< double >::type K(KSEXP);
    Rcpp::traits::input_parameter< double >::type D(DSEXP);
    Rcpp::traits::input_parameter< double >::type N(NSEXP);
    rcpp_result_gen = Rcpp::wrap(logrhonkCalcVarSel(Elogpi, carray, cmatrix, K, D, N));
    return rcpp_result_gen;
END_RCPP
}
// epsCalcVarSel
arma::cube epsCalcVarSel(double K, double maxNCat, double D, double N, arma::mat prioreps, arma::mat X, arma::mat rnk, arma::vec c);
RcppExport SEXP _VICatMix_epsCalcVarSel(SEXP KSEXP, SEXP maxNCatSEXP, SEXP DSEXP, SEXP NSEXP, SEXP priorepsSEXP, SEXP XSEXP, SEXP rnkSEXP, SEXP cSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type K(KSEXP);
    Rcpp::traits::input_parameter< double >::type maxNCat(maxNCatSEXP);
    Rcpp::traits::input_parameter< double >::type D(DSEXP);
    Rcpp::traits::input_parameter< double >::type N(NSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type prioreps(priorepsSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type rnk(rnkSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type c(cSEXP);
    rcpp_result_gen = Rcpp::wrap(epsCalcVarSel(K, maxNCat, D, N, prioreps, X, rnk, c));
    return rcpp_result_gen;
END_RCPP
}
// CpostdeltaCalc
arma::vec CpostdeltaCalc(arma::vec c, double a, double D);
RcppExport SEXP _VICatMix_CpostdeltaCalc(SEXP cSEXP, SEXP aSEXP, SEXP DSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type c(cSEXP);
    Rcpp::traits::input_parameter< double >::type a(aSEXP);
    Rcpp::traits::input_parameter< double >::type D(DSEXP);
    rcpp_result_gen = Rcpp::wrap(CpostdeltaCalc(c, a, D));
    return rcpp_result_gen;
END_RCPP
}
// logeta1Calc
arma::vec logeta1Calc(arma::cube Elogphi, arma::mat rnk, arma::vec Elogdelta, double K, double D, double N);
RcppExport SEXP _VICatMix_logeta1Calc(SEXP ElogphiSEXP, SEXP rnkSEXP, SEXP ElogdeltaSEXP, SEXP KSEXP, SEXP DSEXP, SEXP NSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::cube >::type Elogphi(ElogphiSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type rnk(rnkSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type Elogdelta(ElogdeltaSEXP);
    Rcpp::traits::input_parameter< double >::type K(KSEXP);
    Rcpp::traits::input_parameter< double >::type D(DSEXP);
    Rcpp::traits::input_parameter< double >::type N(NSEXP);
    rcpp_result_gen = Rcpp::wrap(logeta1Calc(Elogphi, rnk, Elogdelta, K, D, N));
    return rcpp_result_gen;
END_RCPP
}
// logeta2Calc
arma::vec logeta2Calc(arma::cube lognullphi, arma::mat rnk, arma::vec Elogminusdelta, double K, double D, double N);
RcppExport SEXP _VICatMix_logeta2Calc(SEXP lognullphiSEXP, SEXP rnkSEXP, SEXP ElogminusdeltaSEXP, SEXP KSEXP, SEXP DSEXP, SEXP NSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::cube >::type lognullphi(lognullphiSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type rnk(rnkSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type Elogminusdelta(ElogminusdeltaSEXP);
    Rcpp::traits::input_parameter< double >::type K(KSEXP);
    Rcpp::traits::input_parameter< double >::type D(DSEXP);
    Rcpp::traits::input_parameter< double >::type N(NSEXP);
    rcpp_result_gen = Rcpp::wrap(logeta2Calc(lognullphi, rnk, Elogminusdelta, K, D, N));
    return rcpp_result_gen;
END_RCPP
}
// cCalc
arma::vec cCalc(arma::vec logeta1, arma::vec clse, double D);
RcppExport SEXP _VICatMix_cCalc(SEXP logeta1SEXP, SEXP clseSEXP, SEXP DSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type logeta1(logeta1SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type clse(clseSEXP);
    Rcpp::traits::input_parameter< double >::type D(DSEXP);
    rcpp_result_gen = Rcpp::wrap(cCalc(logeta1, clse, D));
    return rcpp_result_gen;
END_RCPP
}
// sumDElogphiCalcVarSel
arma::mat sumDElogphiCalcVarSel(arma::cube carray, arma::mat cmatrix, double K, double D, double N);
RcppExport SEXP _VICatMix_sumDElogphiCalcVarSel(SEXP carraySEXP, SEXP cmatrixSEXP, SEXP KSEXP, SEXP DSEXP, SEXP NSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::cube >::type carray(carraySEXP);
    Rcpp::traits::input_parameter< arma::mat >::type cmatrix(cmatrixSEXP);
    Rcpp::traits::input_parameter< double >::type K(KSEXP);
    Rcpp::traits::input_parameter< double >::type D(DSEXP);
    Rcpp::traits::input_parameter< double >::type N(NSEXP);
    rcpp_result_gen = Rcpp::wrap(sumDElogphiCalcVarSel(carray, cmatrix, K, D, N));
    return rcpp_result_gen;
END_RCPP
}
// logrhonkCalcProfCat
arma::mat logrhonkCalcProfCat(arma::vec Elogpi, arma::mat Elogtheta, arma::vec y, arma::cube carray, arma::mat cmatrix, double K, double D, double N);
RcppExport SEXP _VICatMix_logrhonkCalcProfCat(SEXP ElogpiSEXP, SEXP ElogthetaSEXP, SEXP ySEXP, SEXP carraySEXP, SEXP cmatrixSEXP, SEXP KSEXP, SEXP DSEXP, SEXP NSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type Elogpi(ElogpiSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Elogtheta(ElogthetaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::cube >::type carray(carraySEXP);
    Rcpp::traits::input_parameter< arma::mat >::type cmatrix(cmatrixSEXP);
    Rcpp::traits::input_parameter< double >::type K(KSEXP);
    Rcpp::traits::input_parameter< double >::type D(DSEXP);
    Rcpp::traits::input_parameter< double >::type N(NSEXP);
    rcpp_result_gen = Rcpp::wrap(logrhonkCalcProfCat(Elogpi, Elogtheta, y, carray, cmatrix, K, D, N));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_VICatMix_rnkCalc", (DL_FUNC) &_VICatMix_rnkCalc, 4},
    {"_VICatMix_ElogphiCalc", (DL_FUNC) &_VICatMix_ElogphiCalc, 6},
    {"_VICatMix_ElogphiLCalc", (DL_FUNC) &_VICatMix_ElogphiLCalc, 4},
    {"_VICatMix_logrhonkCalc", (DL_FUNC) &_VICatMix_logrhonkCalc, 5},
    {"_VICatMix_epsCalc", (DL_FUNC) &_VICatMix_epsCalc, 7},
    {"_VICatMix_firstepsCalc", (DL_FUNC) &_VICatMix_firstepsCalc, 7},
    {"_VICatMix_CpriorepsCalc", (DL_FUNC) &_VICatMix_CpriorepsCalc, 4},
    {"_VICatMix_CpostepsCalc", (DL_FUNC) &_VICatMix_CpostepsCalc, 4},
    {"_VICatMix_sumDElogphiCalc", (DL_FUNC) &_VICatMix_sumDElogphiCalc, 4},
    {"_VICatMix_priorepsminusoneCalc", (DL_FUNC) &_VICatMix_priorepsminusoneCalc, 4},
    {"_VICatMix_epsminusoneCalc", (DL_FUNC) &_VICatMix_epsminusoneCalc, 4},
    {"_VICatMix_epsminuspriorepsCalc", (DL_FUNC) &_VICatMix_epsminuspriorepsCalc, 5},
    {"_VICatMix_cmatrixCalc", (DL_FUNC) &_VICatMix_cmatrixCalc, 5},
    {"_VICatMix_ElogthetaCalcCat", (DL_FUNC) &_VICatMix_ElogthetaCalcCat, 3},
    {"_VICatMix_betaCalc", (DL_FUNC) &_VICatMix_betaCalc, 6},
    {"_VICatMix_CpriorbetaCalc", (DL_FUNC) &_VICatMix_CpriorbetaCalc, 3},
    {"_VICatMix_CpostbetaCalc", (DL_FUNC) &_VICatMix_CpostbetaCalc, 3},
    {"_VICatMix_respthetaCalc", (DL_FUNC) &_VICatMix_respthetaCalc, 5},
    {"_VICatMix_nullphiCalc", (DL_FUNC) &_VICatMix_nullphiCalc, 5},
    {"_VICatMix_firstbetaCalc", (DL_FUNC) &_VICatMix_firstbetaCalc, 6},
    {"_VICatMix_lognullphiCalc", (DL_FUNC) &_VICatMix_lognullphiCalc, 5},
    {"_VICatMix_logrhonkCalcVarSel", (DL_FUNC) &_VICatMix_logrhonkCalcVarSel, 6},
    {"_VICatMix_epsCalcVarSel", (DL_FUNC) &_VICatMix_epsCalcVarSel, 8},
    {"_VICatMix_CpostdeltaCalc", (DL_FUNC) &_VICatMix_CpostdeltaCalc, 3},
    {"_VICatMix_logeta1Calc", (DL_FUNC) &_VICatMix_logeta1Calc, 6},
    {"_VICatMix_logeta2Calc", (DL_FUNC) &_VICatMix_logeta2Calc, 6},
    {"_VICatMix_cCalc", (DL_FUNC) &_VICatMix_cCalc, 3},
    {"_VICatMix_sumDElogphiCalcVarSel", (DL_FUNC) &_VICatMix_sumDElogphiCalcVarSel, 5},
    {"_VICatMix_logrhonkCalcProfCat", (DL_FUNC) &_VICatMix_logrhonkCalcProfCat, 8},
    {NULL, NULL, 0}
};

RcppExport void R_init_VICatMix(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

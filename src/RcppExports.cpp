// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif


RcppExport SEXP _rcpp_module_boot_stan_fit4spatial_joint_logit_mod();
RcppExport SEXP _rcpp_module_boot_stan_fit4spatial_joint_smooth_mod();
RcppExport SEXP _rcpp_module_boot_stan_fit4spatial_joint_smooth_logit_mod();
RcppExport SEXP _rcpp_module_boot_stan_fit4spatial_joint_smooth_unmatched_mod();
RcppExport SEXP _rcpp_module_boot_stan_fit4spatial_mean_smooth_mod();
RcppExport SEXP _rcpp_module_boot_stan_fit4spatial_mean_smooth_unmatched_mod();

static const R_CallMethodDef CallEntries[] = {
    {"_rcpp_module_boot_stan_fit4spatial_joint_logit_mod", (DL_FUNC) &_rcpp_module_boot_stan_fit4spatial_joint_logit_mod, 0},
    {"_rcpp_module_boot_stan_fit4spatial_joint_smooth_mod", (DL_FUNC) &_rcpp_module_boot_stan_fit4spatial_joint_smooth_mod, 0},
    {"_rcpp_module_boot_stan_fit4spatial_joint_smooth_logit_mod", (DL_FUNC) &_rcpp_module_boot_stan_fit4spatial_joint_smooth_logit_mod, 0},
    {"_rcpp_module_boot_stan_fit4spatial_joint_smooth_unmatched_mod", (DL_FUNC) &_rcpp_module_boot_stan_fit4spatial_joint_smooth_unmatched_mod, 0},
    {"_rcpp_module_boot_stan_fit4spatial_mean_smooth_mod", (DL_FUNC) &_rcpp_module_boot_stan_fit4spatial_mean_smooth_mod, 0},
    {"_rcpp_module_boot_stan_fit4spatial_mean_smooth_unmatched_mod", (DL_FUNC) &_rcpp_module_boot_stan_fit4spatial_mean_smooth_unmatched_mod, 0},
    {NULL, NULL, 0}
};

RcppExport void R_init_VSALM(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

#' Mini Pooling algorithm with covariates
#'
#' This function allows you to compute the average number of assays needed per pool using mini pooling algorithm with covariates.
#' @param v Vector of viral load. 
#' @param s Vector of risk score with the same length of viral load.
#' @param K Pool size, default is 5.
#' @param vf_cut Cutoff for individual viral failure, default is 1000. 
#' @param lod Vector of true VL of those undetectable, default is 0.
#' @param method The method used for the average test number calculation, default is "permutation".
#' @param perm_num The number of permutation to be used for the calculation, default is 100.
#' @return 
#' The average number of assays needed per pool and per subject.
#' @keywords Pooling.
#' @export
#' @examples
#' d = Simdata
#' V = d$VL # Viral Load
#' S = d$S # Risk Score
#' mmpa(V, S, K = 3, perm_num = 3)
#' foo; table(foo)





mmpa = function(v, # vector of true VL
                s, # vector of risk score in same length
                K = 5, # pool size
                vf_cut = 1000, # cutoff for individual viral failure
                lod = 0, # vector of true VL of those undetectable
                method = "permutation", 
                perm_num = 100
){
  ##########################
  n = length(v)
  pool.n = (n%/%K)*K
  permindex = sample(rep(1:n, perm_num), perm_num*pool.n, replace = F)
  v0 = v[permindex]
  s0 = s[permindex]
  impafoo = impa(v0, s0, K, vf_cut, lod = 0)
  cat("For pool size of", K, "the average number of assays needed per pool is", 
      mean(impafoo), ", \ni.e. average number of assays per subject is", 
      mean(impafoo)/K, ".")
  impafoo
}
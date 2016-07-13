#' Mini Pooling Test Number
#'
#' This function allows you to compute the average number of assays needed per pool within each provided pooled cluster.
#' @param v Vector of viral load. 
#' @param s Vector of risk score with the same length of viral load.
#' @param K Pool size, default is 5.
#' @param vf_cut Cutoff for individual viral failure, default is 1000. 
#' @param lod Vector of true VL of those undetectable, default is 0.
#' @return 
#' A vector with length equals the number of pools. Each element indicates the estimated number of tests for the corresponding pool.
#' @keywords Pooling, average number of assays.
#' @export
#' @examples
#' d = Simdata
#' V = d$VL # Viral Load
#' S = d$S # Risk Score
#' impa(V, S, K = 5, vf_cut = 1000, lod = 0)

impa = function(
  v, # vector of VL
  s, # vector of risk score in same length
  K = 5, # pool size
  vf_cut = 1000, # cutoff for individual viral failure
  lod = 0 # vector of true VL of those undetectable
){
  n = length(v)
  if (length(s) != n) {
    warning("V and S have different length!")
    n = min(n, length(s))
  }
  num_pool = n %/% K
  foo0 = n - num_pool*K
  if(foo0 !=0) {
    warning("Last (", foo0, ") observation(s) are dropped to make (", num_pool, ") pools!")
    v = v[1:(num_pool*K)]
    s = s[1:(num_pool*K)]
  }
  ### create a matrix of size  K x num.pool.
  ### i.e. each column is a pool of size K
  #v_mat = matrix(v, nrow = K)
  s_mat = matrix(s, nrow = K)
  ### re-order VL based on the rank of S
  ### s.t. the 1st row has the lowest risk score
  ### and the Kth row the highest
  order_mat = apply(s_mat, 2, order)
  foo = c(order_mat) + rep(0:(num_pool-1), rep(K, num_pool)) * K
  v_mat = matrix(v[foo], nrow = K)
  #print(matrix(s_mat[foo], nrow=K))
  #print(v_mat)
  t_mat = apply(v_mat, 2, cumsum)
  #print(t_mat)
  t0_mat = 0
  if(lod > 0){
    v0_mat = v_mat
    v0_mat[v0_mat >= lod] = 0
    t0_mat = apply(v0_mat, 2, function(x) sum(x)-cumsum(x))
  }
  return(1+apply((t_mat+t0_mat)[-1,]>vf_cut, 2, sum))
}



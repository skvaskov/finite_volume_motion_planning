%  original code
%  multivdtmom.m   Date: 1/25/2016
%  This program computes E[prod_{i=1}^nX_i^{k_i}|a_i<X_i<b_i], where (X_1,..,X_n) 
%  follows a multivariate normal distribution with mean mu and covariance
%  matrix S, and lower and upper truncation limits of a_i and b_i.
%  Input:
%  kappa: power of X_i's
%  a: lower truncation limits
%  b: upper truncation limits
%  mu: E[X]
%  S: Var[X]
%  Output
%  M: E[\prod_{i=1}^nX_i^{I_i}|a_i<X_i<b_i], where I_i ranges from 0 to k_i

%  created by Sean Vaskov on 11/5/2019 to remove normalization factor
function M = multivdtmom_unormalized(kappa,a,b,mu,S)
if size(kappa,1)>1
   kappa = kappa';
end
if size(a,2)>1
   a = a';
end
if size(b,2)>1
   b = b';
end
if size(mu,2)>1
   mu = mu';
end
%
%  Compute the moments of doubly truncated multivariate normal
%
M = recintab(kappa,a,b,mu,S);


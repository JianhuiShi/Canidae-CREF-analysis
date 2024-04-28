function [U, S, V, A1, A2] = robustPCA(X)
%%%
% input:
% X: size: m*n
% output:
% A1: low-rank matrix, size: m*n
% A2: sparse matrix, size: m*n
% SVD: A1=U*S*V'
% S: diagonal, size: m*n
% U: orthogonal matrix, size: m*m
% V: orthogonal matrix, size: n*n
% column vector of U and V is eigenvector
%%%

[A1, A2, ~] = inexact_alm_rpca(X);
[U, S, V] = svd(A1, "econ");

end
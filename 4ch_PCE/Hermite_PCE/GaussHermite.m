function [zeros,weight] = GaussHermite(n,sum_weight)
% This function determines the zeros and weights for the
% Gauss-Hermite quadrature of order n>1, in the interval (-INF, +INF).
%% INPUT: 
%   n          - number of nodes and weights we want to calculate
%   sum_weight - (OPTIONAL)is the weighting function integrated 
%          on the entire support. (i.e. sum(weight)=sum_weight)
%          if we use the probability density function sum_weight=1.
%          (DEFAULT sum_weight=1)
%% OUTPUT:
%   zeros  - nodes of the Gauss Hermite quadrature rule
%   weight - Weights associated to the nodes=zeros of the Gauss Hermite
%          quadrature rule
%
%   written by Luca Fenzi - November 19, 2015
%   Contact: fenzi.luca@gmail.com


% setting the default value of sum_weight if it is not specified.
if nargin==1
    sum_weight=1;
end
 
%% Construction of the Jacobi Matrix
% we use the three terms recurrence formula for monic Hermite polynomials:
% J=tridiag(\sqrt(beta_n),alpha_n,\sqrt(beta_n))
% alpha_0=0; beta_n=n 
beta=sqrt((1:n-1));
J=diag(beta,1)+diag(beta,-1);

%% Determing the zeros (zeros) and weights (weight)
    % - the zeros are the eigenvalues of the Jacobi matrix
    % - the weights can be derived from the corresponding eigenvectors.
    %   weight=squared first component of the eigenvector multiplied by the
    %   integral of the weight function on the support

[eigenvectors, sigma]=eig(J);
zeros=diag(sigma);
weight=(sum_weight*eigenvectors(1,:).^2)';   % We choose the first raw of the matrix eigenvectors
                                             % since eigenvectors(:,n) is
                                             % the nth eigenvector
end


function [zeros,weight] = GaussLaguerre(n,sum_weight)
% This function determines the zeros and weights for the
% Gauss-Laguerre quadrature of order n>1, in the interval [0, +INF).
%% INPUT: 
%   n          - number of nodes and weights we want to calculate
%   sum_weight - (OPTIONAL)is the weighting function integrated 
%          on the entire support. (i.e. sum(weight)=sum_weight)
%          if we use the probability density function sum_weight=1.
%          (DEFAULT sum_weight=1)
%% OUTPUT:
%   zeros  - nodes of the Gauss Laguerre quadrature rule
%   weight - Weights associated to the nodes=zeros of the Gauss Laguerre
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
% J=tridiag(\sqrt(b_n),a_n,\sqrt(b_n))
a=2*(1:n)-1;
b=(1:n-1);
J=diag(a)+diag(b,1)+diag(b,-1);

%% Determing the zeros (zeros) and weights (weight)
    % - the zeros are the eigenvalues of the Jacobi matrix
    % - the weights can be derived from the corresponding eigenvectors.
    %   weight=squared first component of the eigenvector multiplied by the
    %   integral of the weight function on the support
    
[eigenvectors, eigenvalues]=eig(J);
[zeros, sorting]=sort(diag(eigenvalues));  % The eigenvalues are sorted
eigenvectors=eigenvectors(:,sorting)';     % eigenvectors are sorted using 
                                           % the sorting of the eigenvalues
weight=sum_weight*eigenvectors(:,1).^2;             

end


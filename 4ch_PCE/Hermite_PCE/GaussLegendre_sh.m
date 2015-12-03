function [zeros, weight] = GaussLegendre_sh(n)
% This function determines the zeros and weights for the
% Gauss shifted Legendre quadrature of order n>1, in the interval [0, 1].
%% INPUT: 
%   n          - number of nodes and weights we want to calculate
%% OUTPUT:
%   zeros  - nodes of the Gauss shifted Legendre quadrature rule
%   weight - Weights associated to the nodes=zeros of the Gauss shifted Legendre
%          quadrature rule
%% NB 
% This function merges the function GaussLegendre and
% GaussLegendre_changeinterval for this reason we set the the weighting 
% function integrated on the entire support equal to 1.
%% REMARK 
% It is possible to use the three term formula for the shifted Legendre
% polynomials, however this is more expensive than using the three terms
% recurrence formula of the Legendre polynomials and then translate in the
% interval [0,1]
%
%   written by Luca Fenzi - November 19, 2015
%   Contact: fenzi.luca@gmail.com

%% Construction of the Jacobi Matrix
% we use the three terms recurrence formula for monic Hermite polynomials:
% J=tridiag(\sqrt(b_n),0,\sqrt(b_n))
b=(1:n-1)./sqrt(4*(1:n-1).^2-1);
J=diag(b,1)+diag(b,-1);

%% Determing the zeros (zeros) and weights (weight) in [-1,1]
    % - the zeros are the eigenvalues of the Jacobi matrix
    % - the weights can be derived from the corresponding eigenvectors.
    %   w=squared first component of the eirenvector multiplied by the
    %   integral of the weight function on the support
[eigenvectors, eigenvalues]=eig(J);
[zeros, sorting]=sort(diag(eigenvalues));  % The eigenvalues are sorted
eigenvectors=eigenvectors(:,sorting)';     % eigenvectors are sorted using 
                                           % the sorting of the eigenvalues
weight=2*eigenvectors(:,1).^2; 

%% Scale parameters for the interval [0,1].
rescaling=1/2;

% Calculation of the nodes in [0,1]
zeros=rescaling*zeros+1/2;
% Calculation of the weights in [0,1]
weight=rescaling*weight;

     

end

function [x, w] = GaussLegendre_sh(n)
    %% Gauss LEGENDRE formula for the interval [0,1]
    % For the uniform distribution
    % x nodes
    % w weights

%% Building the Jacobi Matrix
%we are going to use the recurrence formula for monic Hermite polynomials:
b=(1:n-1)./sqrt(4*(1:n-1).^2-1);
J=diag(b,1)+diag(b,-1);

%% Determing the zeros (x) and weights (w) in [-1,1]
    % - the zeros are the eigenvalues of the Jacobi matrix
    % - the weights can be derived from the corresponding eigenvectors.
    %   w=squared first component of the eirenvector multiplied by the
    %   integral of the weight function on the support
    
[eigenvectors, eigenvalues]=eig(J);
[x, sorting]=sort(diag(eigenvalues));             %The eigenvalues are sorted
eigenvectors=eigenvectors(:,sorting)';      %eigenvectors are sorted using the sorting of the eigenvalues
w=2*eigenvectors(:,1).^2; 

    %% Calculation of the weights
weight=1/2;

    %% Calculation of the nodes in [0,1]
x=weight*x+1/2;
    %% Calculation of the weights in [0,1]
w=weight*w;

     

end
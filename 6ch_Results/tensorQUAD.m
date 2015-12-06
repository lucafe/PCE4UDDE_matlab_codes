function [Z,w,N_node] = tensorQUAD(l,D)
%tensorQUAD
%   Computes the tensor grid quadrature abscissae and weights
%   on the hypercube [0,1]^D using the Clenshaw-Curtis rule.
%
%   The command [Z,w,N]=tensorQUAD(D,l) 
%   is used inside the function StochasticCollocation in order to construct
%   the grid points and the weights for the uniform distributions in the
%   case 'solver','tensor'
%
%   Input Parameters:
%   D - Stochastic dimension
%   l - order of the integration rule
% 
%   Output Parameters:
%   x - grid points (Chebyshev extrema)
%   w - weights (Clenshaw-Curtis rule)
%   N_node - Number of nodes 
%
%
% The calculation of the weights for the Clenshaw-Curtis rule was inspired
% by the algortithm proposed in: Jorg Waldvogel, "Fast construction of the 
% Fejér and Clenshaw-Curtis quadrature rules," BIT Numerical 
% Mathematics 43 (1), p. 001-018 (2004).
%
%   written by Luca Fenzi - November 3, 2015
%   Contact: fenzi.luca@gmail.com

N_1d=2^l+1;   % Number of nodes in one dimension
N_node=N_1d^D;    %Total number of nodes

%% CONSTRUCTION OF THE INDEXES FOR TENSORIZATION
% For further details see the Matlab forum 
% http://www.mathworks.com/matlabcentral/answers/249168-d-for-loops-concatenated#answer_196354
indexes=zeros(N_node,D);
row=zeros(1,D);
for s=1:N_node
  indexes(s,:)=row;
    %% Updating
    %the alghorithm is counting how many indexes arrived to P-1 (Maximum
    %degree)
    h=0;
    while (h<=D-2)&&(row(h+1)==N_1d-1)
        h=h+1;
    end
    if h>0
        row(1:h)=zeros(1,h);   %Put zero all the values that were P-1
        row(h+1)=row(h+1)+1;   % add one to the next index.
    end
    row(1)=mod(s,N_1d);
end

%% CALCULATION OF THE NODES 1D
f=@(i) unique((i~=0)*cos(pi*2^(-i)*(0:(2^i)))');  % Chebyshev extrema
x1d=(f(l)+1)/2;                                   % nodes in one dimension

w1d=clencurt(2^l+1);                              % Weights in one dimension
w1d=w1d/sum(w1d);                                 % Renormalization for the interval [0,1]


%% Tensor product 
Z=zeros(N_node,D);  
w=zeros(N_node,1); 
for s=1:N_node
       Z(s,:)=x1d(indexes(s,:)+1);   % grid points
       w(s)=prod(w1d(indexes(s,:)+1));   % weights
end


%% Compute 1D Clenshaw-Curtis weights
% Reference: Jorg Waldvogel, "Fast construction of the 
% Fejér and Clenshaw-Curtis quadrature rules," BIT Numerical 
% Mathematics 43 (1), p. 001-018 (2004).

function w=clencurt(N1)

if N1==1
     w=2;
else
    N=N1-1;  c=zeros(N1,1);
    c(1:2:N1,1)=(2./[1 1-(2:2:N).^2 ])'; 
    f=real(ifft([c(1:N1,:);c(N:-1:2,:)]));
    w=2*([f(1,1); 2*f(2:N,1); f(N1,1)])/2;
end


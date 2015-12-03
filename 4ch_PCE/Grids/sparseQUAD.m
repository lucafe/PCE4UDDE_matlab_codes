function [x,w,N_node]=sparseQUAD(l,D)
%sparseQUAD
%   Computes the sparse grid quadrature abscissae and weights
%   on the hypercube [0,1]^D using the Clenshaw-Curtis rule.
%
%   The command [Z,w,N]=sparseQUAD(D,l) 
%   is used inside the function StochasticCollocation in order to construct
%   the grid points and the weights for the uniform distributions in the
%   case 'solver','sparse'
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
% This algorithm is inspired by the code SPQUAD of Greg von Winckel
% available at the link https://people.sc.fsu.edu/~jburkardt/m_src/spquad/spquad.html
% or at the Matlab Central File Exchange:  http://www.mathworks.com/matlabcentral/fileexchange/19063-sparse-grid-quadrature/content/spquad.m
%
% It is also suggested to read the forum whose the algorithm SPQUAD was
% improved, in order to understand all the details 
% http://compgroups.net/comp.soft-sys.matlab/sparse-grid-quadrature/883112
%
% The calculation of the weights for the Clenshaw-Curtis rule was inspired
% by the algortithm proposed in: Jorg Waldvogel, "Fast construction of the 
% Fejér and Clenshaw-Curtis quadrature rules," BIT Numerical 
% Mathematics 43 (1), p. 001-018 (2004).
% 
%   written by Luca Fenzi - November 3, 2015
%   Contact: fenzi.luca@gmail.com

%We calculate the grid points and weights on the hypercube [-1,1]^D and
%then we translate them on [0,1]^D;
bpt=repmat([-1,1]',1,D);

% Length and midpoint in each dimesion
len=diff(bpt); mpt=mean(bpt);

if l==0 % zeroth order case is just the midpoint rule
    x=mpt; w=prod(len);
elseif l~=0 && D==1 % Special 1D case
    x=mpt+len*cos(pi*2^(-l)*(0:(2^l))')/2;
    w=clencurt(2^(l)+1)*len/2;
else
        
    % node and weight for single point grid
    x0=mpt;    w0=prod(len);

    % Construct higher order grid hierarchically 
    for j=1:l
        [x,w]=sparsegridnd(D,j,mpt,len);
        [~,I]=intersect(x,x0,'rows');
        w(I)=w(I)+w0; x0=x; w0=w;
    end

end

% Renormalization of the weights for the interval [0,1]^D
w=w/sum(w);

% Translation of the grid points on the hypercube [0,1]^D
x=(x+1)/2;

% Calculation of the number of nodes used in the grids
N_node=size(x,1);


%% Multidimensional nodes and difference weights
function [x,w]=sparsegridnd(n,ord,mpt,len)

% Generates the n-dimensional sparse grid of order ord on the 
% hyperrectangle with centroid mpt and side lengths len, based
% on Chebyshev-Gauss-Lobatto points.

% Define 1D grid points
p=@(i) unique((i~=0)*cos(pi*2^(-i)*(0:(2^i)))');

% Get configurations of all possible subgrids
v=genindex(n,ord); vmax=1+v(1,n);

% Compute all orders of one-dimensional quadrature rules
P=uaf(@(j) p(j), (0:vmax-1));
Q=uaf(@(j) diffweight(j), (0:vmax));

% Take the union of all possible subgrids 
m=size(v,1);
xw=uaf(@(k) getpts(P,Q,v(k,:),mpt,len), (1:m)',1);
x=xw(:,1:n); w=xw(:,n+1);

% Kludge to deal with small errors introduced by ndgrid
% need a better way to do this
roundn=@(a,n) round(a*10^n)/10^n;
x=roundn(x,15);

% Get unique points and weights
[x,ii,jj]=unique(x,'rows'); rows=length(ii); cols=length(jj);  

% Do node condesation for combining weights
w=sparse(jj,(1:cols)',ones(cols,1),rows,cols)*w;


%% One dimensional difference weights
function dw=diffweight(ii)

if ii==0        
    dw=2;    
elseif ii==1    
    dw=[1;-2;1]/3;
else
    q=@(i) clencurt(2^(i)+1);
    dw=q(ii); dw(1:2:end)=dw(1:2:end)-q(ii-1);
end

%% Compute the subset grid points and difference weights
function [XW]=getpts(P,Q,v,mpt,len)

n=size(v,2);
[x{1:n}]=ndgrid(P{1+v});    [w{1:n}]=ndgrid(Q{1+v});
d=size(v,2);

X=uaf(@(k) mpt(k)+x{k}(:)*len(k)/2, (1:d),1);
W=prod(uaf(@(k) w{k}(:), (1:d),1),2);
XW=[X,W];

%% Find all possible combinations of bases
function v=genindex(n, L1, head)

if n==1
    v = L1;
else
    v = uaf(@(j) genindex(n-1, L1-j, j),(0:L1)',1);
end

if nargin>=3
    v = [head+zeros(size(v,1),1) v];
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

%% Shorthand array function with uniform output
% Optional third argument converts output to matrix
function result=uaf(arg,dex,varargin)
result=arrayfun(arg,dex,'UniformOutPut',false);

if nargin>2
    result=cell2mat(result);
end

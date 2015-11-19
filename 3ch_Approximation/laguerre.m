function l = laguerre(n,x)
% LAGUERRE: compute the Laguerre polynomials of degree n.
%           x optional values where we want to evaluate the Laguerre
%           polynomial of degree n
% This function builds the monic version of the Laguerre polynomials 
% It should be applied with the GaussLaguerre quadrature rule
%
% Inputs:
%   - n is the order of the Laguerre polynomial (n>=0).
%   - x is (optional) values to be evaluated on the resulting Laguerre
%     polynomial function.
% 
% There are two possible outputs:
% 1. If x is omitted then h is an array with (n+1) elements that contains
%    coefficients of each Laguerre polynomial term.
%    E.g. calling l = laguerre(3)
%    the output is l = [1 -9 +18 -6], i.e. x^3 - 3/5x
% 
% 2. If x is given, then l = L_n(x) and l is the same size of x.
%    E.g., L2(x) = x^2 - 4x + 2
%    calling l = laguerre(2,[0 1 2])
%    the output is l = [2 -1 -2]

% The main idea was taken by File Exchange of MATLAB Central at the link
% http://it.mathworks.com/matlabcentral/fileexchange/27746-hermite-polynomials/content/hermite.m
%
%   written by Luca Fenzi - November 19, 2015
%   Contact: fenzi.luca@gmail.com

% check n
if( n<0 ), error('The order of Laguerre polynomial must be greater than or equal to 0.'); end
% again check n is an integer
if( 0~=n-fix(n) ), error('The order of Laguerre polynomial must be an integer.'); end

% call the hermite recursive function.
l = laguerre_rec(n);



% evaluate the Laguerre polynomial function, given x
if( nargin==2 )
    y = l(end) * ones(size(x));
    j = 1;
    for i=length(l)-1:-1:1
        y = y + l(i) * x.^j;
        j = j+1;
    end
    % restore the shape of y, the same as x
    l = reshape(y,size(x));
end


function l = laguerre_rec(n)
% This is the reccurence construction of Laguerre polynomials, i.e.:
%   L_0(x) = 1
%   L_1(x) = x-1
%   L_[n+1](x) = (x-(2n+1)) L_n(x) - (n^2) L_[n-1](x)


if( 0==n ), l = 1;
elseif( 1==n ), l = [1 -1];
else
    % monic term: x_L[n+1]
    l1 = zeros(1,n+1);
    l1(1:n) = laguerre_rec(n-1);
    % beta_n: x
    l2 = zeros(1,n+1);
    l2(2:end) = (2*n-1)*laguerre_rec(n-1);
    % alpha_n
    l3 = zeros(1,n+1);
    l3(3:end) = ((n-1)^2)*laguerre_rec(n-2);
    % L_n+1
    l = l1 - l2 - l3;
    
end
function p = legendre_sh(n,x)
% shifted LEGENDRE: compute the shifted Legendre polynomials of degree n.
%           x optional values where we want to evaluate the shifted Legendre
%           polynomial of degree n
% This function builds the shifted Legendre polynomials 
% It should be applied with the GaussLegendre quadrature rule
%
% Inputs:
%   - n is the order of the shifted Legendre polynomial (n>=0).
%   - x is (optional) values to be evaluated on the resulting shifted Legendre
%     polynomial function.
% 
% There are two possible outputs:
% 1. If x is omitted then h is an array with (n+1) elements that contains
%    coefficients of each shifted Legendre polynomial term.
%    E.g. calling p = legendre_sh(3)
%    the output is p = [20 -30 12 -1], i.e. 20*x^3 -30*x^2+12*x-1
% 
% 2. If x is given, then p = P^*_n(x) and p is the same size of x.
%    E.g., P^*_2(x) = 6*x^2 -6*x+1
%    calling p = legendre_sh(2,[0 1 2])
%    the output is p = [1 1 13]
%
% The main idea was taken by File Exchange of MATLAB Central at the link
% http://it.mathworks.com/matlabcentral/fileexchange/27746-hermite-polynomials/content/hermite.m
%
%   written by Luca Fenzi - November 19, 2015
%   Contact: fenzi.luca@gmail.com

% check n
if( n<0 ), error('The order of Legendre polynomial must be greater than or equal to 0.'); end

% again check n is an integer
if( 0~=n-fix(n) ), error('The order of Legendre polynomial must be an integer.'); end

% call the hermite recursive function.
p = legendre_rec(n);



% evaluate the Legendre polynomial function, given x
if( nargin==2 )
    y = p(end) * ones(size(x));
    j = 1;
    for i=length(p)-1:-1:1
        y = y + p(i) * x.^j;
        j = j+1;
    end
    
    % restore the shape of y, the same as x
    p = reshape(y,size(x));
end


function p = legendre_rec(n)
% This is the reccurence construction of shifted Legendre polynomials, i.e.:
%   P^*_0(x) = 1
%   P^*_1(x) = 2x-1
%   P^*_[n+1](x) = A*x P^*_n(x) -B*P^*_n(x)- C* P*_[n-1](x)

% 0th case
if( 0==n ), p = 1;
% 1st case
elseif( 1==n ), p = [2 -1];
% nth case
else
    % A
    p1 = zeros(1,n+1);
    p1(1:n) = (4*(n-1) +2)/n*legendre_rec(n-1);
    % B
    p2 = zeros(1,n+1);
    p2(2:end) = (2*(n-1) +1)/n*legendre_rec(n-1);
    % C
    p3 = zeros(1,n+1);
    p3(3:end) = ((n-1))/(n)*legendre_rec(n-2);
    % P^*_[n+1]
    p = p1 -p2 - p3;
    
end
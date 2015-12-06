function p = legendre_sh(n,x)
% LEGENDRE: compute the Legendre polynomials of degree n.
%           x optional values where we want to evaluate the Legendre
%           polynomial of degree n
% THIS FUNCTION BUILD THE MONIC VERSION OF LEGENDRE POLYNOMIALS 
% IT SHOULD BE APPLIED WITH GAUSS LEGENDRE quadrature

% Inputs:
%   - n is the order of the Legendre polynomial (n>=0).
%   - x is (optional) values to be evaluated on the resulting Legendre
%     polynomial function.
% 
% There are two possible outputs:
% 1. If x is omitted then h is an array with (n+1) elements that contains
%    coefficients of each Legendre polynomial term.
%    E.g. calling p = legendre(3)
%    will result p = [1 0 -0.6 0], i.e. x^3 - 3/5x
% 
% 2. If x is given, then p = P_n(x) and p is the same size of x.
%    E.g., P2(x) = x^2 - 1/3
%    calling p = legendre(2,[0 1 2])
%    will result p = [-1/3 1/6 3+2/3]

% The main idea was taken by File Exchange of MATLAB Central
% link http://it.mathworks.com/matlabcentral/fileexchange/27746-hermite-polynomials/content/hermite.m

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
% This is the reccurence construction of Legendre polynomials, i.e.:
%   P0(x) = 1
%   P1(x) = x
%   P[n+1](x) = x Pn(x) - (n^2)/(4n^2-1) P[n-1](x)


if( 0==n ), p = 1;
elseif( 1==n ), p = [2 -1];
else
    
    p1 = zeros(1,n+1);
    p1(1:n) = (4*(n-1) +2)/n*legendre_rec(n-1);
    
    p2 = zeros(1,n+1);
    p2(2:end) = (2*(n-1) +1)/n*legendre_rec(n-1);
    
    p3 = zeros(1,n+1);
    p3(3:end) = ((n-1))/(n)*legendre_rec(n-2);
    
    p = p1 -p2 - p3;
    
end
function h = hermite(n,x)
% HERMITE:  compute the probabilistic Hermite polynomials of degree n.
%           x optional values where we want to evaluate the hermite
%           polynomial of degree n
% This function builds the monic version of the probabilistic Hermite polynomials 
% It should be applied with the GaussHermite quadrature rule
%
% Inputs:
%   - n is the order of the probabilistic Hermite polynomial (n>=0).
%   - x is (optional) values to be evaluated on the resulting probabilistic 
%     Hermite polynomial function.
% 
% There are two possible outputs:
% 1. If x is omitted then h is an array with (n+1) elements that contains
%    coefficients of each probabilistic Hermite polynomial term.
%    E.g. calling h = hermite(3)
%    the output is h = [1 0 -3 0], i.e. x^3 - 3x
% 
% 2. If x is given, then h = He_n(x) and h has the same size of x.
%    E.g., He_2(x) = x^2 - 1
%    calling h = hermite(2,[0 1 2])
%    the output is h = [-1 0 3]

% The main idea was taken by File Exchange of MATLAB Central at the link:
% http://it.mathworks.com/matlabcentral/fileexchange/27746-hermite-polynomials/content/hermite.m
%
%   written by Luca Fenzi - November 19, 2015
%   Contact: fenzi.luca@gmail.com

% check n
if( n<0 ), error('The order of Hermite polynomial must be greater than or equal to 0.'); end

% again check n is an integer
if( 0~=n-fix(n) ), error('The order of Hermite polynomial must be an integer.'); end

% call the hermite recursive function.
h = hermite_rec(n);

% evaluate the hermite polynomial function, given x
if( nargin==2 )
    y = h(end) * ones(size(x));
    j = 1;
    for i=length(h)-1:-1:1
        y = y + h(i) * x.^j;
        j = j+1;
    end
    % restore the shape of y, the same as x
    h = reshape(y,size(x));
end


function h = hermite_rec(n)
% This is the reccurence construction of the probabilistic Hermite polynomial, i.e.:
%   He_0(x) = 1
%   He_1(x) = x
%   He_[n+1](x) = x He_n(x) - n He_[n-1](x)

% 0th case
if( 0==n ), h = 1;
% 1st case
elseif( 1==n ), h = [1 0];
% n>1
else
    % He_[n-1](x)
    h1 = zeros(1,n+1);
    h1(1:n) = hermite_rec(n-1);
    % He_[n-2](x)
    h2 = zeros(1,n+1);
    h2(3:end) = (n-1)*hermite_rec(n-2);
    % He_[n](x)
    h = h1 - h2;
    
end

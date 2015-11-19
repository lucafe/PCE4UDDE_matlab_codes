function [zeros, weight] = GaussLegendre_changeinterval(zeros,weight,a,b)
% This Function permits to translate the nodes and weights of the Gauss 
% Legendre quadrature rule evaluated in the interval [-1,1] to an arbitrary
% interval [a,b]
% Change the interval of the nodes and weights of the Lagendre polynomials
%% INPUT: 
%       zeros  - zeros of the Legendre polynomials in [-1,1]
%       weight - weight of the Legendre polynomials in [-1,1]
%       [a,b]  - (OPTIONAL) interval where we want to translate the weights 
%              and the zeros
%              DEFAULT [0,1]
%
%% OUTPUT: 
%       zeros  - zeros of the Legendre polynomials in [a,b]
%       weight - weight of the Legendre polynomials in [a,b]
%
%   written by Luca Fenzi - November 19, 2015
%   Contact: fenzi.luca@gmail.com

% setting the default value of the interval [a,b]=[0,1]
if nargin==1
    a=0;
    b=1;
end

%% RESCALING
rescaling=(b-a)/2;
zeros=rescaling*zeros+(b+a)/2;
weight=rescaling*weight;
end


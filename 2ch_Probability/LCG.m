function [x] = LCG(n)
%% LCG (Linear Congruential Generator)
% is the most famous pseudo random generator
% Input      n - number of pseudo random number that we want to generate
% Output     x - vector containing the sequence of random numbers

% It is initialized with the following parameters:
m=2^(31)-1;      %modulus
c=0;             %increment
a=7^5;           %multiplier


x=zeros(1,n); % initialization of the vector that are going to contain 
              % the pseudo random number sequences
x(1)=randi([1 m-1],1,1); % the seed is chosen randomly

for i=2:n
    x(i)=mod(a*x(i-1)+c,m);  % step of the LCG
end

x=x/m; % this renormalization allow the pseudo random number to be in the interval [0,1]
end


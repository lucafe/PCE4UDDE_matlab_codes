function [indexes,N] = Indexes(P)
% This function creates a vector of indexes for the degrees of the different 
% orthogonal polynomials.
% Every polynomials of the basis is defined using a row of indexes 
% Every value in the row indicate the degree of the correspective
% polynomial.
%% For example indexes(5,:)=[3,0,1] 
% indicate that the 5th polynomial of the basis is calculating multiplying
% the 5th orthogonal polynomial of the 1st parameter with the 1th orthogonal polynomial of the 3rd parameter

%% INPUT:
%      P-1  maximum value inside indexes (i.e maximum degree for a specific
%           direction)
%      D    Stochastic dimension of the problem
D=2;

%% OUTPUT:
%      indexes matrix that for each row indicates the different degree of
%              the polynomials that constitute the basis
%      N       Number of rows that we have in indexes 
%              For the Collocation methods this correspond to the number of
%              nodes in the sparse grid

N       = nchoosek(P-1+D,D); 
indexes = zeros(N,D);
v       = ones(1,D);
for s = 1:N   
  indexes(s, 1:D-1) = -diff(v);
  indexes(s, D)     = v(D) - 1;
    %% Update index vector:
    for k = D:-1:1
      v(k) = v(k) + 1;
      if k > 1
        if v(k) <= v(k - 1);
          break;  % Exit for k loop
        end
        v(k) = 1;
      end
    end
end


end


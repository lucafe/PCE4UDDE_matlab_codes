function [Mean,Variance,N,Prob,Skew,Kurt] = StochasticCollocation(myUDDE,l,D,n_par,varargin)
%StochasticCollocation 
%  Computes specific statistics of a UDDE throught the Stochastic
%  Collocation Method or the Monte Carlo Method.
%  In particular it evaluates the Stochastic collocation using the Lagrange
%  interpolation on sparse grid and tensor grid.
%  the UDDE must be described in a script, see the example: myUDDE.m 
%  In order to construct your own file.
% 
% INPUT 
%  'myUDDE' -    script file where the UDDE is described.
%             This is a Script and it must satisfy the requirement for
%             eigAM.m file and should be in the same folder of this file
%             Moreover it must contain the Chaos Legendre expansion of the random
%             inputs
%
%     l  -    order of the integration rule
%
%     D  -    Stochastic dimension of the problem
%
%     n_par  -   Numbers of parameter of the problem (par>=D) since the number
%             of random parameter is D but it is possible to have a problem
%             with fixed random inputs also. 

% OUTPUT 
%     Mean     -  Mean of the rightmost eigenvalue of the IG matrix
%     Variance -  Variance of the rightmost eigenvalue of the IG matrix
%     N        -  Number of nodes used
%     Prob     -  Probability that the rightmost eigenvalue is positive (ie
%               instable DDE)
%     Skew     -  Skewness of the rightmost eigenvalue (3rd standardized moment) 
%     Kurt     -  Kurtosis of the rightmost eigenvalue (4th standardized moment) 
%
% OPTIONS
%  'eigs'     Number of eigenvalues nearest to 0 calculated using the function eigAM
%             M should be sufficienly large in such a way that the set of
%             the M eigenvalues nearest 0 contains the rightmost
%             eigenvalues of the discrete Infinitesimal Generator      
%             'eigs',20   (DEFAULT)
%
%    'solver'   Selection of the solver that we want to use
%               'sparse'  Lagrange interpolation using smolyak sparse grid
%               'tensor'  Lagrange interpolation using tensor grid
%               'MC'      Monte Carlo method
%               'solver','sparse' (DEFAULT)
%
% To run the function StochasticCollocation must be in the same folder also
% the functions:
% eigAM.m 
% myUDDE.m (or your file with the UDDE that you want to study)
% sparseQUAD.m
% tensorQUAD.m
%
% The file eigAM  is used to solve the random eigenvalue problem at the
% collocation points, it is a deterministic code available in:
% D. Breda, S. Maset and R. Vermiglio, "Stability of linear delay
%      differential equations - A numerical approach with MATLAB", in
%      Control, Automation and Robotics, T. Basar, A. Bicchi and M. Krstic
%      eds., Springer, New York, 2015, ISBN 978-1-4939-2106-5.
%
%   written by Luca Fenzi - November 3, 2015
%   Contact: fenzi.luca@gmail.com




%% DEFINITION OF THE OPTIONS 'eig' and 'solver'
p = inputParser;

% Default Options
defaultEIG=20;
defaultsolver='sparse';

%Possible solver options
expectedsolver={'sparse','tensor','MC'};

%Definition of the field 'eig' and 'grid'
addOptional(p,'eigs',defaultEIG,@isnumeric);
addOptional(p,'solver',defaultsolver,...
    @(x) any(validatestring(x,expectedsolver)));
parse(p,varargin{:});
%p.Results.eigs       is the result of the option 'eigs'
%p.Results.solver     is the result of the option 'solver'

%% CONSTRUCTION OF THE GRID
% Z        vector indicating all the points of the grid
%          Z(i,j) is the jth coordinate of the ith Collocation Point
%          the points are situated in the hypercube [0,1]^D
% w        Weights associated at each collocation points 
%          the weights are requested only using the solver 'sparse' and
%          'tensor'
% N        number of total nodes 

% SPARSE STOCHASTIC COLLOCATION
if strcmp(p.Results.solver,'sparse')==1
    [Z,w,N]=sparseQUAD(D,l);
% TENSOR STOCHASTIC COLLOCATION
elseif strcmp(p.Results.solver,'tensor')==1
    [Z,w,N] = tensorQUAD(l,D);
% MONTE CARLO METHOD
% the grid points are randomly chosen in [0,1]^D
% It is not need to compute the weights
elseif strcmp(p.Results.solver,'MC')==1
    N=(2^l+1)^D; % Number of points
    Z=rand(N,D);
end

%% CONSTRUCTION of the GRID FOR THE RANDOM PARAMETER
%  We translate the initial grid points in the hypercube request from the
%  model myUDDE.m
%  For this reason we are going to extract the PCE of the inputs from the
%  script myUDDE.m

% UDDE MODEL CONVERSION
% Extract the chaos coefficients of the random inputs
[Input_PCE,d]=conversion(myUDDE,n_par);

%Research of fixed parameters
Parameters=Z;
s=0; %Counting the stochastic parameters
for i=1:n_par
   %If there is a fixed parameter the alghorithm adds a column with the
   %fixed value
   if length(Input_PCE{i})==1 
       Parameters=[Parameters(:,1:i), Input_PCE{i}*ones(N,1), Parameters(:,1:end)];
   elseif length(Input_PCE{i})>1 
       %Through the polynomial chaos we reconstruct the sampled grid of the respective parameter.
        s=s+1; %Updating the stochastic parameters
        
        %% RECONSTRUCTION OF THE DISTRIBUTION VIA LEGENDRE POLINOMIALS
        % We use the three term recursion to evaluate this points
        step2=ones(size(Z(:,s)));
        step1=2*Z(:,s)-ones(size(Z(:,s)));
        Parameters(:,i)=Input_PCE{i}(1)*step2+Input_PCE{i}(2)*step1;
            for j=2:length(Input_PCE{i})-1
                a=(4*(j-1) +2)/j*Z(:,s)-(2*(j-1) +1)/j;
                b=((j-1))/(j);
                new=a.*step1-b*step2;
                step2=step1;
                step1=new;
                
                Parameters(:,i)=Input_PCE{i}(j+1)*step1+Parameters(:,i);
            end
   end
end
   
%% EVALUATING THE SOLUTION in the COLLOCATION POINTS
lambda_grid=zeros(N,1);
ro=p.Results.eigs/d;
for h=1:N  
    %Evaluation of the M eigenvalues nearest to 0
    [lambda]=eigAM(myUDDE,Parameters(h,:),p.Results.eigs);
    indice=find(abs(lambda)<ro);
    lambda_grid(h)=max(real(lambda(indice)));   
end
%% POST PROCESSOR
% Evaluation of the stastics from lambda_grid using the different solvers


%% MONTE CARLO
if strcmp(p.Results.solver,'MC')==1 
    %Monte carlo integration
    Mean=(1/N)*sum(lambda_grid);
    Variance=(1/N)*sum((lambda_grid-Mean).^2);
    Prob=(1/N)*sum(lambda_grid>0);
    Skew=sum((lambda_grid-Mean).^3)/(N*Variance^(3/2));
    Kurt=sum((lambda_grid-Mean).^4)/(N*Variance^(2));
% The computation for the statistics of 'solver', 'sparse' and 'tensor' is the same
% Since we already compute the weights for the different grids
else % LAGRANGE STOCHASTIC COLLOCATION
    % Lagrange interpolation
    Mean=sum(lambda_grid.*w);
    Variance=sum( ((lambda_grid-Mean).^2).*w);
    Prob=sum((lambda_grid>0).*w);
    Skew=sum( ((lambda_grid-Mean).^3).*w)/(Variance^(3/2));
    Kurt=sum( ((lambda_grid-Mean).^4).*w)/(Variance^(2));
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% OTHER FUNCTION REQUIRED FOR RUNNING THE MAIN ALGORITHM 
%%  [Mean,Variance,N,Prob,Skew,Kurt] = StochasticCollocation(myUDDE,l,D,n_par)
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% TO CONVERT THE SCRIPT 
function [Input_PCE,d]=conversion(myUDDE,n_par)
par=ones(1,n_par);   % this is just to compile the file
Input_PCE = cell(1, n_par); 
eval(myUDDE);


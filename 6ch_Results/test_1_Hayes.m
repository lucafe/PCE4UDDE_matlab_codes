%% Uncertain Hayes Equation
% Example of a UDDE with stochastic dimension 2
% The probability that the UDDE is unstable is 50% 
% Evaluated using the Hayes stability chart

%% EQUATION
% x'(t)=a*x(t)+b*x(t-1)
% Hence

%% MEMO LIST OF POSSIBLE PARAMETERS
% Atilde(t)=par(1);   
%               Uniform in [-1.5,-0.5]
%               Input_PCE{1}=[-1, 0.5]
% Btilde_1=par(2);
%               Uniform in [0.5,1.5]
%               Input_PCE{2}=[1, 0.5]



%% DIMENSION OF THE DDE
d=1; % it's not a system of equations
     % it's a single equation with two offset parameters
%% STOCHASTIC DIMENSION
D=2;
%% Number of parameters
n_par=2;
     

%% CURRENT TIME TERM
Atilde=@(t,d,par) par(1); 
% Legendre Polynomial chaos expansion of the parameter a.
Input_PCE{1}=[-1, 0.5];

%% DISCRETE DELAY TERMS
dd=1; %INPUT discrete delays row vector dd=[d_{1},...,d_{q}]>=0


Btilde{1}=@(t,d,par) par(2); %INPUT: dxd matrix or call to Btilde1.m
% Legendre Polynomial chaos expansion of the parameter b.
Input_PCE{2}=[1, 0.5];

%% DISTRIBUTED DELAY TERMS
% EMPTY
l=[]; %INPUT left integration extrema row vector l=[l_{1},....,l_{w}]>=0
r=[]; %INPUT right integration extrema row vector r=[r_{1},....,r_{w}]>=0
Ctilde{1}=@(t,theta,d,par) []; %INPUT: dxd matrix or call to Ctilde1.m



%% Uncertain Cushing Hayes Equation
% Example of a UDDE with stochastic dimension 3
% equation (5.6) of the Unpublished article of Rossana Vermiglio

%% EQUATION
%x'(t)=a*x(t)+b*x(t-1)+c*int_-1^0 x(t+theta)exp(5*theta)cos(4*pi*theta) diff theta
%Hence

%% MEMO LIST OF POSSIBLE PARAMETERS
% Atilde(t)=par(1);   
%               Uniform in [-1,0]
%               Input_PCE{1}=[-0.5, 0.5]  
%Btilde_1=par(2);
%               Uniform in [-0.5*(1+exp(1)),0.5*(1-exp(1))]
%               Input_PCE{2}=[-0.5*exp(1), 0.5]  
%Ctilde_1=par(3);
%               Uniform in [-1,1]
%               Input_PCE{3}=[0, 1]  


%d=1            %Dimension of the UDDE

%% DIMENSION OF THE DDE
d=1; % it's not a system of equations
     % it's a single equation with two offset parameters
%% STOCHASTIC DIMENSION
D=3;
%% Number of parameters
n_par=3;
     


%% CURRENT TIME TERM
Atilde=@(t,d,par) par(1); %INPUT: dxd matrix or call to Atilde.m
%Legendre Polynomial chaos expansion of the parameter a.
Input_PCE{1}=[-0.5, 0.5];

%% DISCRETE DELAY TERMS
dd=1; %INPUT discrete delays row vector dd=[d_{1},...,d_{q}]>=0
Btilde{1}=@(t,d,par) par(2); %INPUT: dxd matrix or call to Btilde1.m
Input_PCE{2}=[-0.5*exp(1), 0.5]; %Legendre Polynomial chaos expansion of the parameter b.


%% DISTRIBUTED DELAY TERMS
% EMPTY
l=1; %INPUT left integration extrema row vector l=[l_{1},....,l_{w}]>=0
r=0; %INPUT right integration extrema row vector r=[r_{1},....,r_{w}]>=0
Ctilde{1}=@(t,theta,d,par) exp(5*theta)*par(3)*cos(4*pi*theta); 
% Distributed term
Input_PCE{3}=[0, 1];
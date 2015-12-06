%% INTRAVENOUS GLUCOSE TOLERANCE TEST
%
% This script defines functions and parameters of the systems of 2 linear
% Uncertain Delay Differential Equations (UDDEs) 
% Linearization of the system (4) and (5) in 
% A. De Gaetano and O. Arino, "Mathematical modelling of the intravenous 
% glucose and tolerance test", Journal of Mathematical Biology (2000), VOL
% 40, pages 136 - 168
%
% We linearized the dynamics of glucose G and insuline I in the non trivial
% equilibrium (G;I) and we obtain the following linear UDDE
%
% X' =+xi{1}* X    +xi{2}* Y
% Y' =             +xi{3}* Y       +xi{4}int_{-xi{5}}^0x(t+s)ds
%
% where xi{1},...,xi{5} are random parameters approximated via PCE and 
% evaluated by the script  'data_IGTT.m'
% The PCE of this parameter is saved in the file 
% 'PCE_parameters_IGTT'
% That we recall
load('PCE_parameters_IGTT')
%
% The script can be used in
% -"eigAM.m" to approximate the eigenvalues of the infinitesimal generator
%  in the AUTONOMOUS case for specific parameters.
% -"StochasticCollocation" to approximate the Polynomail Expansion of the
%  rightmost eigenvalue of the infinitesimal generator in the AUTONOMOUS
%  case.
%
% Version 1.1, October 22, 2015. 
%
% References:
%  [1] D. Breda, S. Maset and R. Vermiglio, "Stability of linear delay
%      differential equations - A numerical approach with MATLAB", in
%      Control, Automation and Robotics, T. Basar, A. Bicchi and M. Krstic
%      eds., Springer, New York, 2015, ISBN 978-1-4939-2106-5.
%  [2] D. Breda, S. Maset and R. Vermiglio, "Approximation of eigenvalues
%      of evolution operators for linear retarded functional differential
%      equations", SIAM J. Numer. Anal. 50(3):1456-1483, 2012.
%  [3] D. Breda, S. Maset and R. Vermiglio, "Pseudospectral differencing
%      methods for characteristic roots of delay differential equations",
%      SIAM J. Sci. Comput. 27(2):482-495, 2005.

%% MEMO LIST OF POSSIBLE PARAMETERS
% Atilde=[par(1), par(2); 0, par(3)];
%               Input_PCE{1}=xi{1}
%               Input_PCE{2}=xi{2}
%               Input_PCE{3}=xi{3}
% Btilde_1=par(2);
%               Uniform 
%               Input_PCE{2}=[-0.5*exp(1), 0.5]



%% DIMENSION OF THE DDE
d=2; % it's not a system of equations
     % it's a single equation with two offset parameters
%% STOCHASTIC DIMENSION
D=5;
%% Number of parameters
n_par=5;


%% CURRENT TIME TERM
Atilde=@(t,d,par) [par(1), par(2); 0, par(3)]; %INPUT: dxd matrix or call to Atilde.m
Input_PCE{1}=xi{1};      % Polynomial chaos expansion of the first coefficients
Input_PCE{2}=xi{2};      % Polynomial chaos expansion of the second coefficients
Input_PCE{3}=xi{3};      % Polynomial chaos expansion of the fourth coefficients


%% DISCRETE DELAY TERMS
dd=[]; %INPUT discrete delays row vector dd=[d_{1},...,d_{q}]>=0
Btilde{1}=@(t,d,par) []; %INPUT: dxd matrix or call to Btilde1.m
%Input_PCE{d*d+1}=[];   % Polynomial chaos expansion of the first
                        %coeffiecients of B(1)
%...
%Btilde{q}=@(t,d,par) []; %INPUT: dxd matrix or call to Btildeq.m
%Input_PCE{d*d*(q-1)+1}=[];   % Polynomial chaos expansion of the first 
                              %  coeffiecients of B(q)
%...

%% DISTRIBUTED DELAY TERMS
l=par(5); %INPUT left integration extrema row vector l=[l_{1},....,l_{w}]>=0
Input_PCE{5}=xi{5};
r=0; %INPUT right integration extrema row vector r=[r_{1},....,r_{w}]>=0
Ctilde{1}=@(t,theta,d,par) [0, 0;0, par(4)]; %INPUT: dxd matrix or call to Ctilde1.m
Input_PCE{4}=xi{4}; % Polynomial chaos expansion 

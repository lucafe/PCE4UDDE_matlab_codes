%% Hermite PCE of standard normal distribution 
% (TRIVIAL EXAMPLE)

%% GERM = standard normal:
mu=0;    % Mean of the normal distribution
sigma=1; % Variance of the normal distribution

%Number of nodes used in the Gaussian Quadrature
N=1000;
[x,w]=GaussLegendre_sh(N); % to approximate evaluate the numerator
[xh, wh]=GaussHermite(N);  % to approximate the norm 
                           % (it is possible to evaluate it anaytically)

% The invers ecdf of the standard normal is approximated by
% norminv(p,mu,sigma)                           
                           


%% Calculation of the coefficients of Hermite PCE
P=10; % higher order of the coefficient
xi=zeros(1,P+1);
for i=1:P+1
    norm=sum(wh.*hermite(i-1,xh).*hermite(i-1,xh));
    % Integrand
    f=@(x) norminv(x,mu,sigma).*hermite(i-1,norminv(x,mu,sigma));
    %Gaussian quadrature
    xi(i)=(sum(w.*norminv(x,mu,sigma).*hermite(i-1,norminv(x,mu,sigma))))/norm;
end


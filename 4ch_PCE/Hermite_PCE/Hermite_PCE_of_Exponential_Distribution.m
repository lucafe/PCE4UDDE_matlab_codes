%% Hermite PCE of Exponential distribution 

%% GERM = standard normal:
mu=0;    % Mean of the normal distribution
sigma=1; % Variance of the normal distribution

%Number of nodes used in the Gaussian Quadrature
N=1000;
[x,w]=GaussLegendre_sh(N);   % to approximate evaluate the numerator
[xh, wh]=GaussHermite(N);    % to approximate the norm 
                             % (it is possible to evaluate it anaytically)

% Inverse of the exponential distribution with mean 1:
F=@(x) -log(x);  
% The invers ecdf of the standard normal is approximated by
% norminv(p,mu,sigma)

%% Calculation of the coefficients of Hermite PCE
P=10; % higher order of the coefficient
xi=zeros(1,P+1);
for i=1:P+1
    norm=sum(wh.*hermite(i-1,xh).*hermite(i-1,xh)); 
    Integrand=@(x) F(x).*hermite(i-1,norminv(x,mu,sigma));
    % Gaussian quadrature
    xi(i)=(sum(w.*Integrand(x)))/norm;
end

%% PLOT of the Probability density functions for the different order of PCE
%  Approximation
subplot(1,2,1)
hold on
grid on
% Sample of the germ
n_sample=100000;                 
normal_sample=randn(n_sample,1); 
% range=-1:0.01:10;

PCE=@(x) xi(1);
for i=2:5
    PCE=@(x) xi(i).*hermite(i-1,x)+PCE(x);
    
    ensamble_output=PCE(normal_sample); 
    
    [f,range] = ksdensity(ensamble_output);
    plot(range,f)  %Plot the PDF of i-th order of the distribution obtained by the PCE
end
% Exact PDF of the exponential distribution
plot(range,exppdf(range,1));  
axis([-1, 5, 0, 1])
legend('1st-order','2nd-order','3rd-order','4th-order','exact','Location','southoutside','Orientation','horizontal')
%% SAVE FIGURE
% figname='4PCE_PCE_expXiu_pdf.tex';
% matlab2tikz(figname,'width','1\textwidth') 


%% PLOT of the Chaos Coefficient
subplot(1,2,2)
grid on
plot(0:6,xi(1:7),'-s')
axis([0 6 -1 1])
xlabel('Index')
ylabel('Coefficients')


 set(0,'defaultlinelinewidth',1.5)
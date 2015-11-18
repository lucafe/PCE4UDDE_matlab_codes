%% SAMPLING METHOD FOR THE UNIFORM DISTRIBUTION
% This figure shows the exact pdf compared with the pdf obtained by the
% kernel density estimator through a sample of the uniform distribution
% using the Linear Congruential Genera

% Sample of the uniform distribution
n=1000000;
sample_unif=LCG(n);

%  Parameters of the uniform distribution
a=0;
b=1;

% Support
support=-1:0.01:2;

% Kernel density estimator
[pdf,range] = ksdensity(sample_unif,support);

%  exact pdf
p=@(x) 1/(b-a)*(x>=a & x<=b);

hold on
plot(support,p(support),'-.','LineWidth',1.5);
plot(range,pdf,'-','LineWidth',1.5)
axis([-1 2 -0.2 1.2])
legend('pdf','sample',-1)
title('Sampling a uniform distribution')




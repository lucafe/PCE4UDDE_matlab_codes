%% SAMPLING METHOD FOR THE STANDARD NORMAL DISTRIBUTION
% This figure shows the exact pdf compared with the pdf obtained by the
% kernel density estimator through a sample of the normal distribution
% evaluated by the BOX-MULLER ALGORITHM

set(0,'defaultlinelinewidth',1.2);

% Support
support=-5:0.01:5;

%% Box Muller algorithm
n=10000;
sample_unif=LCG(n)*2*pi;                 % sample of the uniform in [0,2*pi)
sample_exp=-2*log(LCG(n));               % sample of the exponential E(1/2)
r_1=sqrt(sample_exp).*cos(sample_unif);  % First sample of the normal
r_2=sqrt(sample_exp).*sin(sample_unif);  % Second sample of the normal

hold on
plot(support,normpdf(support,0,1),'-.','LineWidth',2) % exact pdf

% Estimated pdf through the first sample
[pdf,range] = ksdensity(r_1,support);
plot(range,pdf,'-','LineWidth',1.5)

% Estimated pdf through the second sample
[pdf,range] = ksdensity(r_2,support);
plot(range,pdf,'-','LineWidth',1.5)

legend('pdf','Y_1','Y_2')

title('Sampling normal distributions')




%% NORMAL DISTRIBUTION N(mu,sigma)
% figures with the pdfs and the cdfs

% Support of the normal distribution
support=-5:0.01:5;

%% Probability density function
% normpdf -> function used to evaluate the exact pdf
hold on
set(0,'defaultlinelinewidth',1.5);
plot(support,normpdf(support,-2,0.2))
plot(support,normpdf(support,0,1))
plot(support,normpdf(support,0,2))
plot(support,normpdf(support,-2,0.5))
legend('\mu=-2, \sigma=0.2', '\mu=0, \sigma=1', '\mu=0, \sigma=2','\mu=-2, \sigma=0.5')
title('Probability density function')


%% Cumulative distribution function
% normcdf -> function used to evaluate the exact cdf
figure
hold on
set(0,'defaultlinelinewidth',1.5);
plot(support,normcdf(support,-2,0.2))
plot(support,normcdf(support,0,1))
plot(support,normcdf(support,0,2))
plot(support,normcdf(support,-2,0.5))
title('Cumulative distribution function')

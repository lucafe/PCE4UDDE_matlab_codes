%% EXPONENTIAL DISTRIBUTION E(lambda)
% figures with the pdfs and the cdfs

% Support of the exponential distribution
support=0:0.01:10;

%% Probability density function
% exppdf -> function used to evaluate the exact pdf
hold on
set(0,'defaultlinelinewidth',1.5);  %width of the lines
plot(support,exppdf(support,0.5))
plot(support,exppdf(support,1))
plot(support,exppdf(support,2))
legend('\lambda=2', '\lambda=1', '\lambda=0.5')
title('Probability density function')


%% Cumulative distribution function
% expcdf -> function used to evaluate the exact cdf
figure
hold on
set(0,'defaultlinelinewidth',1.5);
plot(support,expcdf(support,0.5))
plot(support,expcdf(support,1))
plot(support,expcdf(support,2))
title('Cumulative distribution function')


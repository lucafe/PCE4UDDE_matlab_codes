%% SAMPLING METHOD FOR THE EXPONENTIAL DISTRIBUTION
% This figure shows the exact pdf compared with the pdf obtained by the
% kernel density estimator through a sample of the exponential distribution
%
%   written by Luca Fenzi - November 18, 2015
%   Contact: fenzi.luca@gmail.com

% Since we want that the exact pdf and the pdf obtained using the kernel 
% density estimator have the same color, we require that the default colors
% of Matlab restart cyclically after 3 times.
cmp=colormap(lines(3));
set(groot,'defaultAxesColorOrder',cmp)
% in this way the colors turn ciclically

% Support
support=0:0.01:10;

% Rate of the different exponential distribution
mu=zeros(1,3);
mu(1)=0.5;
mu(2)=1;
mu(3)=2;

% Plotting the exact pdf
hold on
plot(support,exppdf(support,mu(3)),'-.','LineWidth',1.5)
plot(support,exppdf(support,mu(2)),'-.','LineWidth',1.5)
plot(support,exppdf(support,mu(1)),'-.','LineWidth',1.5)
legend('\lambda=0.5', '\lambda=1', '\lambda=2')

% Samples
n=1000000; % number of sample
%% mu(1)
sample_exp=-log(rand(1,n))/mu(1);  % application of the inversion method
[pdf,range] = ksdensity(sample_exp,support); % estimate of the density through the sample 
plot(range,pdf,'-','LineWidth',1.5)


%% mu(2)
sample_exp=-log(rand(1,n))/mu(2);
[pdf,range] = ksdensity(sample_exp,support);
plot(range,pdf,'-','LineWidth',1.5)

%% mu(3)
sample_exp=-log(rand(1,n))/mu(3);
[pdf,range] = ksdensity(sample_exp,support);
plot(range,pdf,'-','LineWidth',1.5)

axis([0 7 0 2])
title('Sampling exponential distributions')

% delete the impostation for the colors
set(groot,'defaultAxesColorOrder','remove')
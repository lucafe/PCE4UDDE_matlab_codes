%% DENSITY ESTIMATORS
% This is an example to ilustrate how it is possible to plot the sampling
% distribution
% In particular the first figure analyse the relationship between the exact
% pdf and that one estrimated through the histogram density estimator
% The second figure compares the exact pdf with the pdfs generated by the
% histogram density estimator and the normal kernel density estimator.
%
%   written by Luca Fenzi - November 18, 2015
%   Contact: fenzi.luca@gmail.com

%% 1st Figure
rng(4)            % Set seed, it is usefull to obtain precise figure 
                  % every time that you run the script
sample=randn(10000,1); % sample from a strandard normal distribution
[hist,range]=histcounts(sample,'Normalization','pdf'); % Histogram density estimator
                  % Specifying 'Normalization','pdf' we require to obtain a
                  % continuos pdf.
centerbin=(range(1:length(range)-1)+range(2:length(range)))/2;  % Recenter the bin if we use 
                  % a piecewise function instead of a histogram to plot the
                  % histogram
hold on
histogram(sample,'Normalization','pdf')          % The histogram 
plot(centerbin,hist,'-r','LineWidth',2)     % The piecewise function 
                                            % obtained by the histogram
plot(range,normpdf(range,0,1),'-.black','LineWidth',1.5) % exact pdf
axis([-4 4 0 0.45])
title('Histogram estimator')

%% 2nd Figure
figure
hold on
[kpdf,rangek] = ksdensity(sample);       % Kernel Density estimator
                                         % calculation
plot(range,normpdf(range,0,1),'-.black','LineWidth',1.5)  % exact pdf
plot(centerbin,hist,'-r','LineWidth',2)                   % pdf obtained by the histogram estimator
plot(rangek,kpdf,'-b','LineWidth',2)                      % pdf obtained by the Kernel density estimator
axis([-4 4 0 0.45])
title('Density estimators')

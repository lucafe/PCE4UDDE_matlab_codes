%% UNIFORM DISTRIBUTION U(a,b)
% figures with the pdf and cdf of a generic unifor distribution
%
%   written by Luca Fenzi - November 18, 2015
%   Contact: fenzi.luca@gmail.com

% Definition of the interval of the uniform distribution 
a=0;
b=1;

% Support
support=-1:0.001:2;

%% Probability density function
p=@(x) 1/(b-a)*(x>=a & x<=b);  % exact pdf

set(0,'defaultlinelinewidth',1.5);
plot(support,p(support))
axis([-1 2 -0.2 1.2])
title('Probability density function')
ax.XTickLabel={' ', ' ',  'a', ' ', 'b', ' ', ' '};
ax.YTickLabel={' ','min = 0', ' ', ' ', ' ', ' ','max = 1/(b-a)'};
ax.XTickLabelRotation = 45;
ax.YTickLabelRotation = 45;

%% Cumulative distribution function
F=@(x) (x-a)./(b-a).*(x>=a & x<=b)+1*(x>b);  % exact cdf
figure
set(0,'defaultlinelinewidth',1.5);
plot(support,F(support))
axis([-1 2 -0.2 1.2])
ax.XTickLabel=[' ', ' ',  'a', ' ', 'b', ' ', ' '];
ax.YTickLabel=[' ','min = 0', ' ', ' ', ' ', ' ','max = 1'];
ax.XTickLabelRotation=45;
ax.YTickLabelRotation=45;
title('Cumulative distribution function')

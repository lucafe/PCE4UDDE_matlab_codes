%% INVERSION METHOD
% This script plot an image that can be used to explain the inversion
% method
%
%   written by Luca Fenzi - November 18, 2015
%   Contact: fenzi.luca@gmail.com

% we use the inverse of the CDF of the  normal distribution
% x=norminv(p,mu,sigma); % inverse of the cdf of the normal
% p=normcdf(x,mu,sigma); % cdf of the normal


mu=0; % mean
sigma=0.5; % standard deviation
support=-3:0.01:3;  %support of the normal


hold on 
set(0,'defaultlinelinewidth',1.2);

plot(support,normcdf(support,mu,sigma),'black','LineWidth',1.5) 
% plot(support,normpdf(support,mu,sigma),'cyan')
%       If you want it is possible to see also the pdf of the distribution
%       even though in this case it is only misleading
% legend('CDF','PDF')
axis([-3 3 -0.1 1.1])

u_v=0.375;   % point U
x_v=norminv(u_v,mu,sigma);  % point X
plot(-3:0.01:x_v,u_v*ones(size(-3:0.01:x_v)),'--blue'); % line from the cdf to X
plot(x_v*ones(size(0:0.01:u_v)),0:0.01:u_v, '--blue');  % line from U to the cdf
% Boundary 
plot(support,zeros(size(support)),'-.red');  % 0 bound
plot(support,ones(size(support)),'-.red');   % 1 bound

text(-3,u_v,'u','FontWeight','bold','HorizontalAlignment','right') % position for the test of U
text(x_v,0,'x','FontWeight','bold','HorizontalAlignment','center') % position for the test of X

title('Inversion Method')




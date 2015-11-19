%% POCHHAMMER SYMBOL figure
% This script plots the first 6 pochhammer symbol function in the domain
% [-4,4]. The pochhammer functions are evaluated using the built in Matlab
% function pochhammer.
%
%   written by Luca Fenzi - November 19, 2015
%   Contact: fenzi.luca@gmail.com

% Domain of the functions
domain=-4:0.1:4;

% PLOT of the Pochhammer symbols
hold on 
for i=0:5
    plot(domain,pochhammer(domain,i),'LineWidth',1.5);
end
axis([-4 4 -7.5 10])
grid on
title('Pochhammer Symbol');
xlabel('$x$','Interpreter','Latex')
ylabel('$(x)_n$','Interpreter','Latex')


%% TO SAVE THE FIGURE IN TIKZ
% figname='figure_pochhammer.tex';
% matlab2tikz(figname,'width','1\textwidth') 
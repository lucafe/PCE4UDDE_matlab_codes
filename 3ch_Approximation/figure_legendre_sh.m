%% SHIFTED LEGENDRE POLYNOMIALS figure
% This script plots the first 5 shifted Legendre polynomials, highlithing the zeros 
% This script use the function legendre_sh to constuct the polynomials
% and the function GaussLegendre_sh to evaluate the zeros
%
%   written by Luca Fenzi - November 19, 2015
%   Contact: fenzi.luca@gmail.com

% Domain of the polynomials
domain=0:0.02:1;

% Since we want that the exact pdf and the pdf obtained using the kernel 
% density estimator have the same color, we require that the default colors
% of Matlab restart cyclically after 5 times.
cmp=colormap(lines(5));
set(groot,'defaultAxesColorOrder',cmp)
% in this way the colors turn ciclically

set(0,'defaultlinelinewidth',1.5);

% PLOT of the shifted Legendre polynomials
hold on
for i=0:4
    plot(domain,legendre_sh(i,domain))
end

% PLOT THE ZEROS
for i=0:4
    x=GaussLegendre_sh(i);
    plot(x,zeros(size(x)),'o')
end

axis([0, 1, -1.2, 1.2])
grid on
y=ylabel('$P_n(x)$','Interpreter','Latex');
l=legend('$P_0(x)$', '$P_1(x)$', '$P_2(x)$', '$P_3(x)$', '$P_4(x)$',...
    'Location', 'southoutside', 'Orientation','horizontal');
set(l,'Interpreter','latex')
t=title('Legendre polynomials');

%remove the color impostation given
set(groot,'defaultAxesColorOrder','remove')

%% TO SAVE THE FIGURE IN TIKZ
% figname='figure_legendrePsh.tex';
% matlab2tikz(figname,'width','1\textwidth') 
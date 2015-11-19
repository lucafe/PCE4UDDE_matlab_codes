%% PROBABILISTIC HERMITE POLYNOMIALS figure
% This script plots the first 5 probabilistic Hermite polynomials, 
% highlithing the zeros.
% This script use the function hermite to constuct the polynomials
% and the function GaussHermite to evaluate the zeros
%
%   written by Luca Fenzi - November 19, 2015
%   Contact: fenzi.luca@gmail.com

% Domain of the polynomials
domain=-2:0.1:2;

% Since we want that the exact pdf and the pdf obtained using the kernel 
% density estimator have the same color, we require that the default colors
% of Matlab restart cyclically after 5 times.
cmp=colormap(lines(5));
set(groot,'defaultAxesColorOrder',cmp)
% in this way the colors turn ciclically

set(0,'defaultlinelinewidth',1.5);

% PLOT of the Hermite polynomials
hold on
for i=0:4
    plot(domain,hermite(i,domain))
end

% PLOT THE ZEROS
for i=0:4
    x=GaussHermite(i);
    plot(x,zeros(size(x)),'o')
end

axis([-2, 2, -8, 4])
grid on
y=ylabel('$He_n(x)$','Interpreter','Latex');
l=legend('$He_0(x)$', '$He_1(x)$', '$He_2(x)$', '$He_3(x)$', '$He_4(x)$',...
    'Location', 'southoutside','Orientation','horizontal');
set(l,'Interpreter','latex')
t=title('probabilistic Hermite polynomials');

%remove the color impostation given
set(groot,'defaultAxesColorOrder','remove')

%% TO SAVE THE FIGURE IN TIKZ
% figname='figure_hermiteHe.tex';
% matlab2tikz(figname,'width','1\textwidth') 
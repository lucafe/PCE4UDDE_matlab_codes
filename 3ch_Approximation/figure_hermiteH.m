%% PHYSICISTS HERMITE POLYNOMIALS figure
% This script plots the first 5 physicists Hermite polynomials.
% This script uses the Matlab function hermiteH to constuct the polynomials
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
    plot(domain,hermiteH(i,domain))
end

axis([-2, 2, -30, 20])
grid on
y=ylabel('$H_n(x)$','Interpreter','Latex');
l=legend('$H_0(x)$', '$H_1(x)$', '$H_2(x)$', '$H_3(x)$', '$H_4(x)$',...
     'Location', 'southoutside','Orientation','horizontal');
set(l,'Interpreter','latex')
t=title('physicists Hermite polynomials');


%remove the color impostation given
set(groot,'defaultAxesColorOrder','remove')

%% TO SAVE THE FIGURE IN TIKZ
% figname='figure_hermiteH.tex';
% matlab2tikz(figname,'width','1\textwidth') 
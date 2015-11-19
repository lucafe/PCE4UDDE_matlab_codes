%% LAGUERRE POLYNOMIALS figure
% This script plots the first 5 Laguerre polynomials, highlithing the zeros 
% This script use the function laguerre to constuct the polynomials
% and the function GaussLaguerre to evaluate the zeros
%
%   written by Luca Fenzi - November 19, 2015
%   Contact: fenzi.luca@gmail.com

% Domain of the polynomials
domain=-1:0.1:5;

% Since we want that the exact pdf and the pdf obtained using the kernel 
% density estimator have the same color, we require that the default colors
% of Matlab restart cyclically after 5 times.
cmp=colormap(lines(5));
set(groot,'defaultAxesColorOrder',cmp)
% in this way the colors turn ciclically

set(0,'defaultlinelinewidth',1.5);

% PLOT of the Laguerre polynomials
hold on
for i=0:4
    plot(domain,laguerre(i,domain))
end


plot(domain,laguerre(0,domain)) % since the 0th laguerre polynomials does not
                                % have any zero we replot the function in
                                % such a way that the colors start again
                                % ciclycally.
                              
% PLOT THE ZEROS
for i=1:4
    x=GaussLaguerre(i);
    plot(x,zeros(size(x)),'o')
end

axis([0, 5, -20, 15])
grid on
y=ylabel('$L_n(x)$','Interpreter','Latex');
l=legend('$L_0(x)$', '$L_1(x)$', '$L_2(x)$', '$L_3(x)$', '$L_4(x)$',...
    'Location', 'southoutside', 'Orientation','horizontal');
set(l,'Interpreter','latex')
t=title('Laguerre polynomials');

%remove the color impostation given
set(groot,'defaultAxesColorOrder','remove')

%% TO SAVE THE FIGURE IN TIKZ
% figname='figure_laguerreL.tex';
% matlab2tikz(figname,'width','1\textwidth') 
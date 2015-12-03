%% PCE of the solution of the stochastic equation:
% y'(t)=xi*y(t),
% y(0)=1,
% where xi is a standard normal distribution

% The solution of the stochastic ODE is evaluated using the Galerkin
% projection
% Even though it is possible to evaluate Mean and Variance analytically:
%% Exact Mean 
Ex_Mean=@(x) exp(x.^2/2);
%% Exact Variance
Ex_Var=@(x) exp(x.^2).*(exp(x.^2)-1);

%% OFFSET OF THE CODE
Max_expansion=6; % maximum order of the Hermite PCE +1
N=1000;  %Gauss-Hermite quadrature approximation
[x,w]=GaussHermite(N);
%% Inner product with 3 Hermite polynomials
H3=zeros(Max_expansion,Max_expansion);
for i=1:Max_expansion
    for j=1:Max_expansion                
        H3(i,j)=sum(w.*hermite(1,x).*hermite(j-1,x).*hermite(i-1,x));
    end
end

%% Normalization coefficient of the hermite polynomials
H2=zeros(1,Max_expansion);
for i=1:Max_expansion
    H2(i)=factorial(i-1);                  
end

%% Initialization of the errors 
err_mean=zeros(1,Max_expansion);
err_var=zeros(1,Max_expansion);                 

%% CALCULATION OF THE CHAOS COEFFICIENTS
for P=2:Max_expansion;
    save('coefficients.mat', 'H2', 'H3', 'P')
    
    %% Initial value
    Initial_y=zeros(1,P);
    Initial_y(1)=1;
    %% Resolution of the Model using ode23
    [T,Y] = ode23(@stochastic_ode_PCE,[0 1],Initial_y);
    
    %% Calculation of the errors
    err_mean(P-1)=max(abs(Ex_Mean(T)-Y(:,1)));
    var_approx=(Y(:,2:P).^2)*H2(2:P)';
    err_var(P-1)=max(abs(Ex_Var(T)-var_approx));
end
%% FIGURES

% size of my screen obtained using maxwinsize=get(gcf,'position');
% maxwinsize=[1 45 1280 689]; 
% set(gcf,'position',maxwinsize);
%% 1st FIGURE = Convergence PLOT
figure
subplot(1,2,1)
%hold on
semilogy(1:Max_expansion,err_mean,'-o');
hold on
semilogy(1:Max_expansion,err_var,'-o');
legend('Mean','Variance','Location','southoutside','Orientation','horizontal'); %legend('Mean','Variance');
xlabel('N-1')
title('Error')
% figname='4PCE_SODE_error.tex';
% matlab2tikz(figname,'width','1\textwidth')

%% 2nd Figure = Chaos COEFFICIENTS
hold on 
for i=1:Max_expansion
    plot(T,Y(:,i))
end
plot(T,Ex_Var(T),':');
plot(T,Ex_Mean(T),':');
legend( 'y0', 'y1', 'y2', 'y3', 'y4', 'Exact-variance','Exact-mean','Location','southoutside','Orientation','horizontal')
xlabel('Time')
title('Chaos Coefficients')
% figname='4PCE_SODE_PCE.tex';
% matlab2tikz(figname,'width','1\textwidth')

%% 3rd Figure = Graphics stochastic process
% This graphics plot the different distibution depending on t
% So we use the plot3 function
figure
% set(gcf,'position',maxwinsize);  % Settings for LATEX
% set(gca, 'FontName', 'Arno Pro')
% set(gca,'fontsize', 18);
hold on
grid on
n_sample=1000;
normal_sample=randn(n_sample,1); %Sample the input
Number_Bins=600;

for t=2:length(T)
    PCE=@(x) Y(t,1);   %PCE for every t
    for i=2:Max_expansion
        PCE=@(x) Y(t,i).*hermite(i-1,x)+PCE(x);
    end
    
    ensamble_output=PCE(normal_sample);
    [pdf,range] = ksdensity(ensamble_output);
    plot3(T(t)*ones(1,length(pdf)),range,pdf,'- r')  %Plot the PDF of ith order


end
plot3(T,Ex_Mean(T),zeros(length(T),1),'b');
view(37.5, 30)   % sets the default three-dimensional view
xlabel('Time')
ylabel('Support')
axis([0 1 -2 4 0 10])
title('Pdf of the solution')
% figname='4PCE_SODE_pdf.tex';
% matlab2tikz(figname,'width','1\textwidth')

%% 4th Figure =  Graphic 3D Stochastic Collation
% This function use PCE(x,t)
% and then plot the 3d graphics with surf
figure
% set(gcf,'position',maxwinsize);  % Settings for LATEX
% set(gca, 'FontName', 'Arno Pro')
% set(gca,'fontsize', 18);
hold on
grid on
n_sample=1000;
normal_sample=randn(n_sample,1); %Sample the input

%PCE depending on x in Support
%                 t in [0,1] contained in T
PCE=@(x,t) Y(t,1);
for i=2:Max_expansion
    PCE=@(x,t) Y(t,i).*hermite(i-1,x)+PCE(x,t);
end
range=-2:0.1:4;
Surface=zeros(length(range),length(T));
for i=1:length(T)
    ensamble_output=PCE(normal_sample,i); 
    [Surface(:,i),range]=ksdensity(ensamble_output,range);
end
hold on 
surf(T(1:length(T)),range,Surface,'FaceColor','interp',...
   'EdgeColor','none',...
   'FaceLighting','gouraud')
view(37.5, 30)   % sets the default three-dimensional view
xlabel('Time')
ylabel('Support')
title('Pdf of the solution')
%colorbar
axis([0 1 -2 4 0 10])
PDF=[0,10];
camlight right
colormap(hsv)
caxis(PDF)


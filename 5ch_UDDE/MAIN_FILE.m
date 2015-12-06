%% Global variable used in UDDE_logistic to define the systems of DDEs
load('coefficients.mat')

Max_expansion=6; %max degree for the PCE
                 %Defined in Coefficients_HErmite

                 

%%Initial value
Initial_x=zeros(1,N);
Initial_x(1)=1;
delay=ones(1,N);
Time_interval=[0 1.5];
ddeset('MaxStep', {0.01*abs(Time_interval(1)-Time_interval(2))},...
    'RelTol',1e-4);
sol = dde23(@Hayes_UDDE,delay,Initial_x,Time_interval);
                    %t int  %initial condition

T=sol.x;
X=sol.y;
                    
%% FIGURE
figure
%maxwinsize=[1 45 1280 689]; %size of my screen obtained using maxwinsize=get(gcf,'position');
%set(gcf,'position',maxwinsize);
set(0,'DefaultAxesFontSize',12)
set(0,'defaultlinelinewidth',1.5)


%% Chaos COEFFICIENTS
plot(T,X)
l=legend('$x_0(t)$', '$x_1(t)$', '$x_2(t)$', '$x_3(t)$', '$x_4(t)$','$x_5(t)$',...
    'Location','southoutside','Orientation','horizontal' );
set(l,'Interpreter','Latex');
xlabel('Time')
ylabel('Chaos Coefficients')
%title('Chaos Coefficients');
matlab2tikz('6_res_Hayes_PCE.tex','width','\linewidth')

set(gca, 'FontName', 'Arno Pro')

%% SETTINGS FOR 3d GRAPHICS
% POINT OF VIEW
View =[30, 37.5]; %[30, 37.5]
PDF=[0,8];


%% 3d Graphics
% This graphics plot the different distibution depending on t
% So we use the plot3 function


figure
hold on
grid on
n_sample=100000;
uniform_sample_x=rand(n_sample,1); %Sample the input
uniform_sample_y=rand(n_sample,1);
Number_Bins=600;
min_supp=1;
max_supp=1;

for t=2:length(T)
    %We do not evaluate the ensample at t=1 since the kernel density
    %extimation of an ensample with only one values is a normal (but it
    %should be a degenerate distribution
    
    PCE=@(x,y) X(1,t);   %PCE for every t
    for i=2:N
        PCE=@(x,y) X(i,t).*legendre_sh(indexes(i,1),x).*legendre_sh(indexes(i,2),y)+PCE(x,y);
    end
    
    ensamble_output=PCE(uniform_sample_x,uniform_sample_y);   %sampling of the output distribution
    [pdf,range]=ksdensity(ensamble_output);
    plot3(T(t)*ones(1,length(pdf)),range,pdf,'- r')  %Plot the PDF of ith order
    min_supp=min([min(range), min_supp]);
    max_supp=max([max(range), max_supp]);
end
plot3(T,X(1,:),zeros(length(T),1),'b');   % APPROXIMATED MEAN



view(View)   % sets the default three-dimensional view

xlabel('Time')
ylabel('Support')
zlabel('pdf')
Support=[min_supp-1, max_supp+1];
axis([Time_interval Support PDF])
%axis tight
set(gca, 'FontName', 'Arno Pro')
set(gca,'fontsize', 18);

%% 2nd Graphic 3D
% This function use PCE(x,t)
% and then plot the 3d graphics with surf/mesh
figure
hold on
grid on
n_sample=10000;
uniform_sample_x=rand(n_sample,1); %Sample the input zeta1
uniform_sample_y=rand(n_sample,1); %Sample the input zeta2

%PCE depending on x in Support
%                 t in [0,1] contained in T
PCE=@(x,y,t) X(1,t);
for i=2:N
    PCE=@(x,y,t) X(i,t).*legendre_sh(indexes(i,1),x).*legendre_sh(indexes(i,2),x)+PCE(x,y,t);
end

range=Support(1):0.01:Support(2);
Surface=zeros(length(range),length(T)-1);
for i=2:length(T)
    ensamble_output=PCE(uniform_sample_x,uniform_sample_y,i); 
    [Surface(:,i-1),range]=ksdensity(ensamble_output,range);
    
    %[Surface(:,i)]=histcounts(ensamble_output,range,'Normalization','pdf');
end
surf(T(2:length(T)),range,Surface,'FaceColor','interp',...
   'EdgeColor','none',...
   'FaceLighting','gouraud')
view(View)   % sets the default three-dimensional view
camlight left
xlabel('Time')
ylabel('Support')
zlabel('pdf')
colormap(hsv)   %autumn(5)
caxis(PDF)
axis([Time_interval Support PDF])
set(gca, 'FontName', 'Arno Pro')
set(gca,'fontsize', 18);

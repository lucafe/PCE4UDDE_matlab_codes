% Hermite PCE of Y_1=\alpha+exp(\beta Z) 
% with Z normally distributed
%% PARAMETERS OF THE EQUATIONS
alpha=0;
beta=1;
% Initialization
P=10;  % Number of chaos coefficient 

%% EXACT SOLUTION
% Evaluated by Field and Grigoriu
yi=zeros(1,P);
yi(1)=alpha+exp(beta^2/2);
yi(2:P)=exp(beta^2/2)*beta.^(1:P-1)./factorial(1:P-1);
hold on 
grid on
plot(0:P-1,yi,'-sr')

%% APPROXIMATION
%% Initialization    
N=1000; %Number of nodes used for Gaussian quadrature
[xh, wh]=GaussHermite(N,1);  % Calculation of weights and nodes

%% approximated Hermite PCE 
xi=zeros(1,P);
for i=1:P
  norm=sum(wh.*hermite(i-1,xh).*hermite(i-1,xh));  
                   % It is possible to evaluete it using the exact values
                   % factorial(i-1)/(2^(i-1))
  % Evaluation of the Chaos Coefficients
  xi(i)=(sum(wh.*(alpha+exp(xh*beta)).*hermite(i-1,xh)))/(norm);
end
plot(0:P-1,xi,':*b')
text(0.5,1.7,'Y');
xlabel('Index')
ylabel('Coefficients')
legend('Exact', 'Approximated','Location','southoutside','Orientation','horizontal')
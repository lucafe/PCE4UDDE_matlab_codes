% Polynomial Chaos expansion of Y_2= abs(Z) 
% where Z is normal distribution

% Initialization
P=10;  % Number of chaos coefficient 
%% EXACT SOLUTION
% Evaluated by Field and Grigoriu
yi=zeros(1,P);
for i=1:2:P
    j=(i-1)/2;
    s=0;
    for h=1:j
        s=s+(-1)^h*2^(j-2*h)*factorial(j-h)/(factorial(h)*factorial(2*j-2*h));
    end
    yi(i)=(2/sqrt(2*pi))*(2^j*factorial(j)/factorial(2*j)+s);
end
hold on 
grid on
plot(0:P-1,yi,'-or')


%% APPROXIMATION
%% Initialization    
N=1000; %Number of nodes used for Gaussian quadrature
[xh, wh]=GaussHermite(N,1);  % Calculation of weights and nodes

%% approximated Hermite PCE 
xi=zeros(1,P);
for i=1:P
  norm=sum(wh.*(hermite(i-1,xh).^2));  
                   % It is possible to evaluete it using the exact values
                   % factorial(i-1)/(2^(i-1))
  % Evaluation of the Chaos Coefficients
  xi(i)=(sum(wh.*abs(xh).*hermite(i-1,xh)))/(norm);
end
plot(0:P-1,xi,':*b')

xlabel('Index')
ylabel('Coefficients')
legend('Exact', 'Approximated','Location','southoutside','Orientation','horizontal')
text(0.5,0.7,'V');
%% SAVE FIGURE
% figname='4PCE_Grigoriu.tex';
% matlab2tikz(figname,'width','1\textwidth') 
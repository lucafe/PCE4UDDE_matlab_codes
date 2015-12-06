function dxdt=Hayes_UDDE(t,x,Z)

%   This is a simple example of resolution of the logistic UDDE
%   Using the Galerkin method
%
%   The differential equations:
%
%    x'_l(t) = a_0*x_l(t)+b_0x_l(t-1)+
%             + (sum_i a_1*x_j(t) <P_1(x_1)*Psi_i,Psi_l>/Norm_psi(l))
%             + (sum_i b_1*x_j(t-1) <P_1(x_2)*Psi_i,Psi_l>/Norm_psi(l))
%    x'_l(t) = a_0*x_l(t)+b_0x_l(t-1)+
%             + (sum_i a_1*x_i(t) L3x(l,i)/Norm_psi(l))
%             + (sum_i b_1*x_i(t-1) L3y(l,i)/Norm_psi(l))

%% The original UDDE was 
% y'(t)=alpha*y(t)+beta*y(t-1)

% alpha  and beta are random input
%% Definition random input
% min_alpha=-5.5; max_alpha=-4.5; % alpha is U(min_alpha,max_alpha)
% min_beta=0.5; max_beta=1.5;     % beta is U(min_beta,max_beta)

% %Chaos Coefficients of the random inputs
% b_0=min_beta;
% b_1=(max_beta-min_beta);
% a_0=min_alpha;
% a_1=(max_alpha-min_alpha);
a_0=-5;
b_0=1;
a_1=0.5;
b_1=0.5;

%% P, L2, L3 evaluated in Coefficients_Legendre
load('coefficients.mat')


%there is only one delay in the equation but in the system they are N
%So we use Z(:,i) to indicate y_i(t-tau)

for l=1:N
    lag{l}=Z(:,l);
end




%Solution
%dydt=zeros(P,P);
dxdt=zeros(N,1);
for l=1:N
    for i=1:N
        dxdt(l)=a_1*x(i)*L3x(l,i)+b_1*lag{i}(i)*L3y(l,i)+dxdt(l);
    end 
    dxdt(l)=(1/Norm_psi(l))*(dxdt(l))+a_0*x(l)+b_0*lag{l}(l);
end

end

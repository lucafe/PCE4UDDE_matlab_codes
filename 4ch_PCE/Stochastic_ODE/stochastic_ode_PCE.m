function dy = stochastic_ode_PCE(t,y)
% Function readed by ode23 in the main file

%% P, H2, H3 evaluated in Coefficients_Hermite
load('coefficients.mat')

% derivatives
dy=zeros(P,1);
for i=1:P
    for j=1:P
        dy(i)=y(j)*H3(i,j)+dy(i);
    end 
    dy(i)=-(1/H2(i))*(dy(i));
end
end


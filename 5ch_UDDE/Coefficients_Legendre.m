%% This function shuld be run only once and then the results are saved in the file coefficients.mat file

%Check if the program works:
%        iterationString = ['Iteration #',int2str(i+j)];
%        disp(iterationString)

%% MAXIMUM PCE expansion
P=3;
[indexes,N] = Indexes(P);

%% Normalization coefficient of the hermite polynomials
%% NORM FOR the Legendre Polynomials
norm_legendre=@(n) 1./(2*n+1);
Norm_psi=zeros(1,N);
for l=1:N;
    Norm_psi(l)=norm_legendre(indexes(l,1))*norm_legendre(indexes(l,2));
end

zeros_gauss=1000;  %Gauss-Legendre quadrature approximation
[x,w]=GaussLegendre_sh(zeros_gauss);

%% Inner product with 3 Legendre polynomials case x
L3x=zeros(N,N);
for l=1:N
    for i=1:N 
        % CALCULATION FOR ZETA1
        L3x(l,i)=sum(w.*legendre_sh(1,x).*legendre_sh(indexes(i,1),x).*legendre_sh(indexes(l,1),x));
        % CALCULATION FOR ZETA2
        L3x(l,i)=L3x(l,i)*sum(w.*legendre_sh(indexes(i,2),x).*legendre_sh(indexes(l,2),x));
    end
end

%% Inner product with 3 Legendre polynomials case x
L3y=zeros(N,N);
for l=1:N
    for i=1:N 
        % CALCULATION FOR ZETA2
        L3y(l,i)=sum(w.*legendre_sh(1,x).*legendre_sh(indexes(i,2),x).*legendre_sh(indexes(l,2),x));
        % CALCULATION FOR ZETA2
        L3y(l,i)=L3y(l,i)*sum(w.*legendre_sh(indexes(i,1),x).*legendre_sh(indexes(l,1),x));
    end
end

save('coefficients.mat', 'L3x', 'L3y', 'Norm_psi', 'N','indexes')
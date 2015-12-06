function [lambda,M]=eigAM(myDDE,par,M)
%EIGAM infinitesimal generator approach for linear AUTONOMOUS DDEs.
%  [lambda,M]=EIGAM(myDDE,par,M) approximates the eigenvalues lambda of the
%  infinitesimal generator A associated to the system of d linear and 
%  AUTONOMOUS Delay Differential Equations (DDEs) defined in the script
%  "myDDE.m". It uses pseudospectral differentiation methods to discretize
%  the generator as described in [1, Chapter 5] and originally in [2],
%  according to the notation in [1]. par is a vector of possible model
%  parameters, see "myDDE.m"; M is the discretization index, see [1].
%
%  Example of call: [lambda,M]=eigAM('myDDE',[],M);
%
%  Version 1.0, june 5, 2014.
%
%  References:
%  [1] D. Breda, S. Maset and R. Vermiglio, "Stability of linear delay
%      differential equations - A numerical approach with MATLAB", in
%      Control, Automation and Robotics, T. Basar, A. Bicchi and M. Krstic
%      eds., Springer, New York, 2015, ISBN 978-1-4939-2106-5.
%  [2] D. Breda, S. Maset and R. Vermiglio, "Approximation of eigenvalues
%      of evolution operators for linear retarded functional differential
%      equations", SIAM J. Numer. Anal. 50(3):1456-1483, 2012.
%  [3] L.N. Trefethen, "Spectral methods in Matlab", SIAM, 2000.

%% DDE MODEL CONVERSION
[d,tau,p,A,B,C,flagC]=conversion(myDDE,par);
Id=eye(d);

%% PIECEWISE MESH DEFINITION
dtau=diff(tau)/2;           %half-length of delay intervals
Mk=em1(2*dtau,p,M);         %piecewise spectral discretization of delay
                            %intervals: Mk(k)+1 nodes in
                            %[-tau_{k},-tau_{k-1}]
                            %em1: M points in total
                            %em2: M points on the largest interval
MMk=[0,cumsum(Mk)];         %cumulative Mk's
M=MMk(end);                 %M+1=total number of nodes in [-max(tau),0]

%% MATRIX BLOCK-DIAGONAL: DIFFERENTIATION PART
AM=zeros(M+1);                                %initialization
for k=1:p                                     %along [-tau_{k},-tau_{k-1}]
    OmegaMk=dtau(k)*cos((0:Mk(k))*pi/Mk(k))-...
        (tau(k)+tau(k+1))/2;                  %Chebyshev II nodes
    Dk=difmat(OmegaMk);                       %differentiation matrix
    AM(MMk(k)+2:MMk(k+1)+1,...
        MMk(k)+1:MMk(k+1)+1)=Dk(2:Mk(k)+1,:); %update
end                                           %
AM=kron(AM,Id);                               %block structure for systems

%% MATRIX 1ST BLOCK-ROW: SPLICING CONDITION
AM(1:d,1:d)=A(1,d,par);                    %update with current time matrix
for k=1:p                                     %along [-tau_{k},-tau_{k-1}]
    Bk=B{k};                                  %discrete delay matrix
    AM(1:d,MMk(k+1)*d+1:(MMk(k+1)+1)*d)=...
        Bk(1,d,par);                          %update
    if flagC==1                               %distributed delays
        Ck=C{k};                              %distributed delay kernel
        OmegaMk=dtau(k)*cos((0:Mk(k))*pi/Mk(k))-...
            (tau(k)+tau(k+1))/2;              %Chebyshev II nodes
        wqk=quadwei(Mk(k));                   %quadrature weights
        for j=0:Mk(k)                         %quadrature addends
            AM(1:d,(MMk(k)+j)*d+1:(MMk(k)+j+1)*d)=...
                AM(1:d,(MMk(k)+j)*d+1:(MMk(k)+j+1)*d)+dtau(k)*wqk(j+1)*...
                Ck(1,OmegaMk(j+1),d,par);     %update
        end                                   %
    end                                       %
end                                           %

%% COMPUTATION OF EIGENVALUES
lambda=eig(AM);                       %approximated eigenvalues
[~,ind]=sort(real(lambda),'descend'); %ordered by decreasing real part
lambda=lambda(ind);

%% AUXILIARY FUNCTIONS
function [d,tau,p,A,B,C,flagC]=conversion(myDDE,par)
%CONVERSION convert DDE models.
%  [d,tau,p,A,B,C,flagC]=CONVERSION(myDDE) converts the original system of
%  d linear DDEs
%
%  x'(t)=Atilde(t)x(t)
%   +sum_{u=1}^{q}Btilde_{u}(t)*x(t-d_{u})
%   +sum_{v=1}^{w}int_{-l_{v}}^{-r_{v}}Ctilde_{v}(t,theta)x(t+theta)dtheta
%
%  contained in the script "myDDE.m" into
%
%  x'(t)=A(t)x(t)
%   +sum_{k=p}^{p}B_{k}(t)*x(t-tau_{k})
%   +sum_{k=1}^{p}int_{-tau_{k}}^{-tau_{k-1}}C_{k}(t,theta)x(t+theta)dtheta
%
%  with ordered delays according to the notation in [1, Chapter 7].
%
%  Version 1.0, may 28, 2014.

%% DDE FUNCTIONS AND PARAMETERS
eval(myDDE); %load functions and parameters of the original DDE contained
             %in the script "myDDE.m".

%% DELAYS
q=length(dd); %number of original discrete delay terms
w=length(l);  %number of original distributed delay terms
if w>0                           %eliminate distributed terms with
    ind=find(l-r==0);            %null integration interval
    w=w-length(ind);             %
    Cbar=cell(w,1);              %
    for v=1:w                    %
        if v==ind                %
            Cbar{v}=Ctilde{v+1}; %
            l(v)=l(v+1);         %
            r(v)=r(v+1);         %
        else Cbar{v}=Ctilde{v};  %
        end                      %
    end                          % 
    l=l(1:w);                    %
    r=r(1:w);                    %
    Ctilde=Cbar;                 %
end                              %
tau=unique(sort([0,dd,l,r])); %vector of ordered delays
p=length(tau)-1;              %[tau_{0},tau_{1},...,tau_{p}]

%% CURRENT TIME AND DISCRETE DELAY TERMS
Bbar=cell(p+1,1);
Bbar{1}=@(t,d,par) Atilde(t,d,par);
for k=2:p+1
    Bbar{k}=@(t,d,par) zeros(d);
end
for u=1:q
    k=find(tau==dd(u));
    Bbark=Bbar{k};
    Btildeu=Btilde{u};
    Bbar{k}=@(t,d,par)Bbark(t,d,par)+Btildeu(t,d,par);
end
A=Bbar{1};
B=cell(p,1);
for k=1:p
    B{k}=Bbar{k+1};
end

%% DISTRIBUTED DELAY TERMS
C=cell(p,1);
if isempty(l)
    flagC=0;
else flagC=1;
    for k=1:p
        C{k}=@(t,theta,d,par) zeros(d);
    end
    LR=[l;r];
    [LR,ind]=sort(LR,'descend');
    l=LR(1,:);
    r=LR(2,:);
    ind=find(ind(1,:)==2);
    for v=1:length(ind)
        Ctildev=Ctilde{ind(v)};
        Ctilde{ind(v)}=@(t,theta,d,par)-Ctildev(t,theta,d,par);
    end
    for v=1:w
        Ctildev=Ctilde{v};
        ind=find((-l(v)<=-tau)&(-tau<-r(v)));
        for k=1:length(ind)
            Ck=C{ind(k)-1};
            C{ind(k)-1}=@(t,theta,d,par)Ck(t,theta,d,par)+...
                Ctildev(t,theta,d,par);
        end
    end
end

function Mk=em1(d,p,M)
%EM1 piecewise spectral distribution.
%  Mk=EM1(d,p,M) distributes a total of sum(Mk) points over a grid of p
%  intervals of length d(1),...,d(p) proportionally to spectral accuracy,
%  i.e., (d(i)/Mk(i))^Mk(i) is almost constant for all i=1,...,p, under the
%  constraint that Mk(1),...,Mk(p) are the minimum integers s.t.
%  sum(Mk)>=M.
%
%  Version 1.0, may 28, 2014.

if isempty(d)
    Mk=0;
    return
end
Mk=zeros(1,p);
[~,ind]=max(d);
f=@(x,t,m)(t/x)^x-(d(ind)/m)^m;
m=ceil(M/p);
while sum(Mk)<M
    for k=1:p
        mm=m;
        while (f(mm,d(k),m)<0)&&(mm>=1)
            mm=mm-1;
        end
        Mk(k)=mm;
    end
    m=m+1;
end

function Mk=em2(d,p,M)
%EM2 piecewise spectral distribution.
%  Mk=EM2(d,p,M) distributes a total of sum(Mk) points over a grid of p
%  intervals of length d(1),...,d(p) proportionally to spectral accuracy,
%  i.e., (d(i)/Mk(i))^Mk(i) is almost constant for all i=1,...,p, under the
%  constraint max(Mk)=M.
%
%  Version 1.0, may 28, 2014.

if isempty(d)
    Mk=0;
    return
end
Mk=zeros(1,p);
[~,ind]=max(d);
f=@(x,t)(t/x)^x-(d(ind)/M)^M;
for k=1:p
    mm=M;
    while (f(mm,d(k))<0)&&(mm>=1)
        mm=mm-1;
    end
    Mk(k)=mm;
end

function D=difmat(x)
%DIFMAT differentiation matrix on Chebyshev II nodes.
%  D=DIFMAT(x) returns the differentiation matrix D relevant to the N+1
%  Chebyshev II nodes x according to [3].
%
%  Original version from [3].

N=length(x)-1;
if N==0
    D=0;
    return
end
c=[2;ones(N-1,1);2].*(-1).^(0:N)';
X=repmat(x',1,N+1);
dX=X-X';
D=(c*(1./c)')./(dX+(eye(N+1)));
D=D-diag(sum(D'));

function wq=quadwei(N)
%QUADWEI Clenshaw-Curtis quadrature weights.
%  wq=QUADWEI(N) returns the weights wq of the Clenshaw-Curtis quadrature,
%  i.e., the one based on Chebyshew II nodes in [-1,1], see [3].
%
%  Original version from [3].

if N==0
    wq=2;
else p=pi*(0:N)'/N;
    wq=zeros(1,N+1);
    ii=2:N;
    v=ones(N-1,1);
    if mod(N,2)==0
        wq(1)=1/(N^2-1);
        wq(N+1)=wq(1);
        for k=1:N/2-1
            v=v-2*cos(2*k*p(ii))/(4*k^2-1);
        end
        v=v-cos(N*p(ii))/(N^2-1);
    else
        wq(1)=1/N^2;
        wq(N+1)=wq(1);
        for k=1:(N-1)/2
            v=v-2*cos(2*k*p(ii))/(4*k^2-1);
        end
    end
    wq(ii)=2*v/N;
end
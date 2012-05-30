clear all
N=10;

varphi = 2;
sigma=0.01;
%phi =gprnd(1,varphi,0,N,1);
%phi=sort(phi);
%phi =[ 0.0560 0.1440 0.5907 0.8725 1.1852 1.7711 1.8066 8.4151 12.4015 19.3174];
%x = exp(sigma*randn(N,1));
%x=  [0.9897 1.0134 0.9958 0.9986 1.0090 0.9970 1.0103 0.9966 1.0102 1.0063];


%% Grid
phiv=[0.1:0.01:0.3];
xv=[0.9:0.01:1.1];
% 
% phi=;
% x=1.01;


nobs=length(xv)*length(phiv);



%% Parameters

sig=6;
rho=(sig-1)/sig;
w=1;
xi=0.1;

%% Defined the function
p=zeros(length(phiv),length(xv));

for i=1:length(phiv);
    for j=1:length(xv);
        
    phi=phiv(i);
    x=xv(j);
    
    ptar=w./(rho*phi); %target price
    phix=phi.*x; %product of phi and x

    etaprim = @(p) xi*sig*(2/sig * (1./ptar-1./p)-(p./ptar-1).^2) ; %Quadratique cost

    fun = @(p) 1/rho * w ./ phix - etaprim(p)/(sig-1).* p.^(1+sig); %price rule

    %% Solver Parameters
    pinit=1/rho * w ./ phix;%ones(1,nobs)+5
    options=optimset('Display','Final','MaxFunEvals',50000000);

    p(i,j)=fsolve(fun,pinit,options);
    end;
end;
clear all
close all
N=10;

varphi = 2;
sigma=0.01;
%phi =gprnd(1,varphi,0,N,1);
%phi=sort(phi);
%phi =[ 0.0560 0.1440 0.5907 0.8725 1.1852 1.7711 1.8066 8.4151 12.4015 19.3174];
%x = exp(sigma*randn(N,1));
%x=  [0.9897 1.0134 0.9958 0.9986 1.0090 0.9970 1.0103 0.9966 1.0102 1.0063];


%% Grid
phi=[0.1:0.1:3];
x=[0.9:0.01:1.1];
% 
% phi=;
% x=1.01;
nphi=length(phi);
nx=length(x);
nobs=nx*nphi;

xx=x(repmat(1:length(x),1,length(phi)));

ind=1:nphi;
ind=ind(ones(1,nx),:);
phiphi= phi(:,ind);






%% Parameters

sig=6;
rho=(sig-1)/sig;
w=1;
xi=0.1;

%% Defined the function

ptar=w./(rho.*phiphi); %target price
phix=phiphi.*xx; %product of phi and x

etaprim = @(p) xi*sig*(2/sig * (1./ptar-1./p)-(p./ptar-1).^2) ; %Quadratique cost

fun = @(p) 1/rho * w ./ phix - etaprim(p)/(sig-1).* p.^(1+sig)-p; %price rule

%% Solver Parameters
pinit=1/rho * w ./ phix;%ones(1,nobs)+5
options=optimset('Display','Iter');

p=fsolve(fun,pinit,options);
devia = (p - 1/rho * w ./ phix);

p=reshape(p,nx,nphi);
devia=reshape(devia,nx,nphi);

figure(1);
plot(phi,p(11,:),'k')
hold on;
plot(phi,p(1,:),'b')
plot(phi,p(21,:),'r')
legend('No shock','bad shock','good shock')
xlabel('phi')

figure(2);
plot(phi,devia(11,:),'k')
hold on;
plot(phi,devia(1,:),'b')
plot(phi,devia(21,:),'r')
legend('zero','bad shock','good shock')
xlabel('phi')

figure(3)
plot(phi,abs(devia(16:21,:)))
xlabel('phi')
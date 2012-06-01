clear all
close all


print_price = 'yes'; % 'yes/no' If you want to plto the figure regarding price

%% Grid
phi=[0.1:0.05:3];
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
f=2;

%% Defined the function

ptar=w./(rho.*phiphi); %target price
phix=phiphi.*xx; %product of phi and x

%etaprim = @(p) xi*sig*(2/sig * (1./ptar-1./p)-(p./ptar-1).^2).*p.^(-1-sig);            %Quadratique cost * quantity produce
%etaprim = @(p) xi*2* (p./ptar-1).*1./ptar;                                             %Quadratique cost
etaprim = @(p) xi./(phix).*sig.*(2/sig * (1./ptar-1./p)-(p./ptar-1).^2).*p.^(-1-sig); %Quadratique cost * quantity produce /pdty


%FOC
fun = @(p) 1/rho * w ./ phix - etaprim(p)/(sig-1).* p.^(1+sig)-p; %price rule

%% Solver Parameters
pinit=1/rho * w ./ phix;%ones(1,nobs)+5
options=optimset('Display','Iter');

p=fsolve(fun,pinit,options);
devia = (p - 1/rho * w ./ phix);

p_vec=reshape(p,nx,nphi);
devia=reshape(devia,nx,nphi);

switch print_price
    case 'yes'
figure(1);
plot(phi,p_vec(11,:),'k')
hold on;
plot(phi,p_vec(1,:),'b')
plot(phi,p_vec(21,:),'r')
legend('No shock','bad shock','good shock')
xlabel('phi')
title('Price set')

figure(2);
plot(phi,devia(11,:),'k')
hold on;
plot(phi,devia(1,:),'b')
plot(phi,devia(21,:),'r')
legend('zero','bad shock','good shock')
xlabel('phi')
title('Deviasion from flexible price (p-w/(rho*phi*x))')

figure(3);
plot(phi,abs(devia(11,:)),'k')
hold on;
plot(phi,abs(devia(1,:)),'b')
plot(phi,abs(devia(21,:)),'r')
legend('zero','bad shock','good shock')
xlabel('phi')
title('Deviasion from flexible price (p-w/(rho*phi*x))')

figure(4)
plot(phi,abs(devia(16:21,:)))
xlabel('phi')
title('Deviasion from flexible price |p-w/(rho*phi*x)| for different shock levels')

    case 'no'
end;

%% Quantity

Q=1;
P=2;

q = Q .* (p./P).^(-sig);
r= p.*q;

%% Profit

eta = @(p) xi./phix .* q .* (p./ptar - 1).^2;

pro= r - w .*(f + q./phix)-eta(p);

pro_vec=reshape(pro,nx,nphi);

figure(5);
plot(phi,pro_vec(11,:),'k')
hold on;
plot(phi,pro_vec(1,:),'b')
plot(phi,pro_vec(21,:),'r')
legend('zero shock','bad shock','good shock')
plot(phi,0.*phi,'k--')
xlabel('phi')
xlim([0, 1])
title('Profit')


clear all;
close all;

%% Parameters
phi = 2;
alpha=0.6;
w=0.3;

%Variance of ability shocks
sigma=0.01;

%number of firms
%N=100;
N=[10:10: 1000];



%number of simulations
M=1000;

RES=zeros(3,length(N));
for k=1:length(N)
%% Computing statistics of the economy
varphi = gprnd(1,phi,0,N(k),M);
x = exp(sigma*randn(N(k),M));

y_1 = varphi .* (alpha/w)^(alpha/(1-alpha));
y_idea=varphi .* (alpha/w)^(alpha/(1-alpha));

y_2 = varphi.* x .* (alpha/w)^(alpha/(1-alpha));

Y_1 = sum(y_1);
Y_2 = sum(y_2);

GDP= (Y_2 - Y_1)./Y_1;

%disp('Empirical mean of GDP')
Emp_mean=mean(GDP);

%disp('Empirical std of GDP')
Emp_std=std(GDP);
%disp('Theoretical std of GDP')
Theo_std= mean(sqrt(sum( sigma^2 * (y_idea).^2))./(Y_1));


RES(1,k)=Emp_mean;
RES(2,k)=Emp_std;
RES(3,k)=Theo_std;
RES(4,k)=Emp_std - Theo_std;
end;
err_std=RES(2,:);

%beta = lscov(x,y)

beta = lscov([(1./log(N))',ones(length(N),1)],err_std')
fit=[(1./log(N))',ones(length(N),1)]*beta;

%figure(1);
%subplot(211);
plot(N,RES(2,:),'LineWidth',2);
hold on;
plot(N,fit,'-r','LineWidth',2);
xlabel('Number of firms','fontsize',12,'fontweight','b');
ylabel('Empirical standard deviation','fontsize',12,'fontweight','b');
title('Empirical standard deviation','fontsize',12,'fontweight','b');

% subplot(212);
% plot(N,RES(2,:),'LineWidth',2);
% hold on;
% plot(N,RES(3,:),'-k','LineWidth',2);
% xlabel('Number of firms','fontsize',12,'fontweight','b');
% ylabel('Theoretical standard deviation','fontsize',12,'fontweight','b');
% title('Theoretical vs Empirical standard deviation');
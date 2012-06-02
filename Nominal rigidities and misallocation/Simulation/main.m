
%% Set up


close all

addpath('lib');

%% Parameters
w=1;
sig=4;
rho=(sig-1)/sig;


xi_v=[0:0.01:1];

Pxi=zeros(size(xi_v));
Mxi=zeros(size(xi_v));
p_vec=zeros(length(xi_v),length(phi));

disp('----------------with adjustment cost qty/pdty-------------------')
for i=1:length(xi_v);
    
xi=xi_v(i);
disp(xi);

f=1;

fe=10;
delta=0.1;

L=1;

params = [w sig rho xi f fe delta L];


%% Grid for x and phi

smin = [0.01 0.9];
smax = [0.5 1.1];

orders = [20 20];

cdef=fundefn('lin', orders, smin, smax);
nodes = funnode(cdef);
grid = gridmake(nodes);

ns = size(grid,1);


phi=grid(:,1);
x=grid(:,2);

%% Distribution of phi and x
pdfphi = @(phi) gppdf(phi,1,2,0);
pdfx =@(x) normpdf(x,1,0.1);

weight = pdfphi(phi).*pdfx(x);
weight = weight ./sum(weight);

%weight = 1./ns + 0*weight;

%% Solve for price
disp('Solve for price')
    % Solving info
    options=optimset('Display','Final','TolFun',1e-8);
    
    phix=phi.*x;
    pinit= (w./(rho.*phix));
    

    fun = @(p) RES_quadra_qty_pdty(p,grid,params);
    %fun = @(p) RES_quadra_qty(p,grid,params);
    %fun = @(p) RES_quadra_qty_x_pdty(p,grid,params);
    
    p=fsolve(fun,pinit,options);
    
    [coeff,B]=funfitxy(cdef, grid, p);
    p_fun = @(s) funeval(coeff, cdef, s);
    p_vec(i,:) = p;
    
%% Solve for the other equation of the model

if i>1;
    Minit=Mxi(i-1);
    Pinit=Pxi(i-1);
else
    Minit=10;
    Pinit=20;
end;
%Rinit=1;
eqinit= [Minit Pinit];


disp('Solve for the equilibirum')

fun_eq = @(eq)  RES_other_qty_pdty(eq(1),eq(2),grid,p_fun,weight,params); %qty/pdty
%fun_eq = @(eq)  RES_other_qty(eq(1),eq(2),grid,p_fun,weight,params); %qty
%fun_eq = @(eq)  RES_other_qty_x_pdty(eq(1),eq(2),grid,p_fun,weight,params); %qty*pdty

eq=fsolve(fun_eq,eqinit,options) ;


Mxi(i)=eq(1);
Pxi(i)=eq(2);
R=L;


end;

Pxi_norm=Pxi./Pxi(1);
Mxi_norm=Mxi./Mxi(1);
Pindex= Pxi.^(1-sig)./Mxi;


p_opti=(w./(rho.*phix));
devia=p_vec'-p_opti*ones(1,length(xi_v));


figure(1);
subplot(231)
plot(xi_v,1./Pxi_norm);
xlabel('xi');
title('TFP with Adjustement cost with qty/pdty')

subplot(232)
plot(xi_v, Mxi_norm);
xlabel('xi');
title('M (normalized) with Adjustement cost with qty/pdty')

subplot(233)
plot(xi_v, Pindex);
xlabel('xi');
title('P^(1-sig)/M (normalized) with Adjustement cost with qty/pdty')

subplot(234)
plot(phi(1:20), devia((1:20),1),'k');
hold on;
plot(phi(1:20), devia((1:20),length(xi_v)),'r');
plot(phi(1:20), devia((1:20),floor(length(xi_v)/2)),'b');
xlabel('phi');
title(' p-popti with Adjustement cost with qty/pdty and bad shock')
legend('xi=0','xi high','xi median')


subplot(235)
plot(phi(381:400), devia((381:400),1),'k');
hold on;
plot(phi(381:400), devia((381:400),length(xi_v)),'r');
plot(phi(381:400), devia((381:400),floor(length(xi_v)/2)),'b');
xlabel('phi');
title(' p-popti with Adjustement cost with qty/pdty and good shock')
legend('xi=0','xi high','xi median')

disp('----------------with adjustment cost qty-------------------')
for i=1:length(xi_v);
    
xi=xi_v(i);
disp(xi);

f=1;

fe=10;
delta=0.1;

L=1;

params = [w sig rho xi f fe delta L];


%% Grid for x and phi

smin = [0.01 0.9];
smax = [0.5 1.1];

orders = [20 20];

cdef=fundefn('lin', orders, smin, smax);
nodes = funnode(cdef);
grid = gridmake(nodes);

ns = size(grid,1);


phi=grid(:,1);
x=grid(:,2);

%% Distribution of phi and x
pdfphi = @(phi) gppdf(phi,1,2,0);
pdfx =@(x) normpdf(x,1,0.1);

weight = pdfphi(phi).*pdfx(x);
weight = weight ./sum(weight);

%weight = 1./ns + 0*weight;

%% Solve for price
disp('Solve for price')
    % Solving info
    options=optimset('Display','Final','TolFun',1e-8);
    
    phix=phi.*x;
    pinit= (w./(rho.*phix));
    

    %fun = @(p) RES_quadra_qty_pdty(p,grid,params);
    fun = @(p) RES_quadra_qty(p,grid,params);
    %fun = @(p) RES_quadra_qty_x_pdty(p,grid,params);
    
    p=fsolve(fun,pinit,options);
    
    [coeff,B]=funfitxy(cdef, grid, p);
    p_fun = @(s) funeval(coeff, cdef, s);
    p_vec(i,:) = p;
    
%% Solve for the other equation of the model

if i>1;
    Minit=Mxi(i-1);
    Pinit=Pxi(i-1);
else
    Minit=10;
    Pinit=20;
end;
%Rinit=1;
eqinit= [Minit Pinit];


disp('Solve for the equilibirum')

%fun_eq = @(eq)  RES_other_qty_pdty(eq(1),eq(2),grid,p_fun,weight,params); %qty/pdty
fun_eq = @(eq)  RES_other_qty(eq(1),eq(2),grid,p_fun,weight,params); %qty
%fun_eq = @(eq)  RES_other_qty_x_pdty(eq(1),eq(2),grid,p_fun,weight,params); %qty*pdty

eq=fsolve(fun_eq,eqinit,options) ;


Mxi(i)=eq(1);
Pxi(i)=eq(2);
R=L;


end;

Pxi_norm=Pxi./Pxi(1);
Mxi_norm=Mxi./Mxi(1);
Pindex= Pxi.^(1-sig)./Mxi;


p_opti=(w./(rho.*phix));
devia=p_vec'-p_opti*ones(1,length(xi_v));

figure(2);
subplot(231)
plot(xi_v,1./Pxi_norm);
xlabel('xi');
title('TFP with Adjustement cost with qty')

subplot(232)
plot(xi_v, Mxi_norm);
xlabel('xi');
title('M (normalized) with Adjustement cost with qty')

subplot(233)
plot(xi_v, Pindex);
xlabel('xi');
title('P^(1-sig)/M (normalized) with Adjustement cost with qty')

subplot(234)
plot(phi(1:20), devia((1:20),1),'k');
hold on;
plot(phi(1:20), devia((1:20),length(xi_v)),'r');
plot(phi(1:20), devia((1:20),floor(length(xi_v)/2)),'b');
xlabel('phi');
title(' p-popti with Adjustement cost with qty and bad shock')
legend('xi=0','xi high','xi median')


subplot(235)
plot(phi(381:400), devia((381:400),1),'k');
hold on;
plot(phi(381:400), devia((381:400),length(xi_v)),'r');
plot(phi(381:400), devia((381:400),floor(length(xi_v)/2)),'b');
xlabel('phi');
title(' p-popti with Adjustement cost with qty and good shock')
legend('xi=0','xi high','xi median')

disp('----------------with adjustment cost qty*pdty-------------------')
for i=1:length(xi_v);
    
xi=xi_v(i);
disp(xi);

f=1;

fe=10;
delta=0.1;

L=1;

params = [w sig rho xi f fe delta L];


%% Grid for x and phi

smin = [0.01 0.9];
smax = [0.5 1.1];

orders = [20 20];

cdef=fundefn('lin', orders, smin, smax);
nodes = funnode(cdef);
grid = gridmake(nodes);

ns = size(grid,1);


phi=grid(:,1);
x=grid(:,2);

%% Distribution of phi and x
pdfphi = @(phi) gppdf(phi,1,2,0);
pdfx =@(x) normpdf(x,1,0.1);

weight = pdfphi(phi).*pdfx(x);
weight = weight ./sum(weight);

%weight = 1./ns + 0*weight;

%% Solve for price
disp('Solve for price')
    % Solving info
    options=optimset('Display','Final','TolFun',1e-8);
    
    phix=phi.*x;
    pinit= (w./(rho.*phix));
    

    %fun = @(p) RES_quadra_qty_pdty(p,grid,params);
    fun = @(p) RES_quadra_qty(p,grid,params);
    %fun = @(p) RES_quadra_qty_x_pdty(p,grid,params);
    
    p=fsolve(fun,pinit,options);
    
    [coeff,B]=funfitxy(cdef, grid, p);
    p_fun = @(s) funeval(coeff, cdef, s);
    p_vec(i,:) = p;
    
%% Solve for the other equation of the model

if i>1;
    Minit=Mxi(i-1);
    Pinit=Pxi(i-1);
else
    Minit=10;
    Pinit=20;
end;
%Rinit=1;
eqinit= [Minit Pinit];


disp('Solve for the equilibirum')

%fun_eq = @(eq)  RES_other_qty_pdty(eq(1),eq(2),grid,p_fun,weight,params); %qty/pdty
fun_eq = @(eq)  RES_other_qty(eq(1),eq(2),grid,p_fun,weight,params); %qty
%fun_eq = @(eq)  RES_other_qty_x_pdty(eq(1),eq(2),grid,p_fun,weight,params); %qty*pdty

eq=fsolve(fun_eq,eqinit,options) ;


Mxi(i)=eq(1);
Pxi(i)=eq(2);
R=L;


end;

Pxi_norm=Pxi./Pxi(1);
Mxi_norm=Mxi./Mxi(1);
Pindex= Pxi.^(1-sig)./Mxi;


p_opti=(w./(rho.*phix));
devia=p_vec'-p_opti*ones(1,length(xi_v));

figure(3);
subplot(231)
plot(xi_v,1./Pxi_norm);
xlabel('xi');
title('TFP with Adjustement cost with qty*pdty')

subplot(232)
plot(xi_v, Mxi_norm);
xlabel('xi');
title('M (normalized) with Adjustement cost with qty*pdty')

subplot(233)
plot(xi_v, Pindex);
xlabel('xi');
title('P^(1-sig)/M with Adjustement cost with qty*pdty')

subplot(234)
plot(phi(1:20), devia((1:20),1),'k');
hold on;
plot(phi(1:20), devia((1:20),length(xi_v)),'r');
plot(phi(1:20), devia((1:20),floor(length(xi_v)/2)),'b');
xlabel('phi');
title(' p-popti with Adjustement cost with qty*pdty and bad shock')
legend('xi=0','xi high','xi median')


subplot(235)
plot(phi(381:400), devia((381:400),1),'k');
hold on;
plot(phi(381:400), devia((381:400),length(xi_v)),'r');
plot(phi(381:400), devia((381:400),floor(length(xi_v)/2)),'b');
xlabel('phi');
title(' p-popti with Adjustement cost with qty*pdty and good shock')
legend('xi=0','xi high','xi median')
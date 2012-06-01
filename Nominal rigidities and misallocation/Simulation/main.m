
%% Set up


close all

addpath('lib');

%% Parameters
w=1;
sig=4;
rho=(sig-1)/sig;


xi=0.3;

f=1;

fe=10;
delta=0.1;

L=1;

params = [w sig rho xi f fe delta L];


%% Grid for x and phi

smin = [0.01 0.95];
smax = [0.7 1.05];

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
    options=optimset('Display','Iter');
    
    phix=phi.*x;
    pinit= (w./(rho.*phix));
    

    fun = @(p) RES_quadra_qty_pdty(p,grid,params);
    p=fsolve(fun,pinit,options);
    
    [coeff,B]=funfitxy(cdef, grid, p);
    p_fun = @(s) funeval(coeff, cdef, s);
    
%% Solve for the other equation of the model

Minit=10;
Pinit=20;
%Rinit=1;
eqinit= [Minit Pinit];


disp('Solve for the equilibirum')
fun_eq = @(eq)  RES_other(eq(1),eq(2),grid,p_fun,weight,params);
eq=fsolve(fun_eq,eqinit,options) ;


M=eq(1);
P=eq(2);
R=L;
pro=profit(grid,P,R, p_fun, params);




function pro=profit(s,P,R, p_fun, params)
    

    %% Idio component
    phi=s(:,1);
    x=s(:,2);
    
    %% Parameters
    w=params(1);
    sig=params(2);
    rho=params(3);
    xi= params(4);
    f=params(5);
    fe=params(6);
    delta=params(7);
    L=params(8);
    
    %% Some others useful intermediate computation
    phix=phi.*x;
    ptar= (w./(rho.*phi));
    Q=R/P;
   
    

    %% Profit
    pro = Q.*P^sig./(sig-1).*p_fun([phi,x]).^(-sig) .*( w./phix -xi./phix.*(p_fun([phi,x])./ptar-1 ).^2 -2*xi./phix .*(1./ptar-1./p_fun([phi,x])) ) -w.*f;
    
    
    
end
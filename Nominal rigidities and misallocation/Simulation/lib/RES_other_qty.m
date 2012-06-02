function RES = RES_other_qty_x_pdty(M,P,grid,p_fun,weight,params)

    %% Idio component
    phi=grid(:,1);
    x=grid(:,2);
    
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

    %% revenu = size of the economy
    R=L;
    
    %% Solve for mu
    
    mu = (profit_qty(grid,P,R, p_fun, params)>0).*weight ./ sum((profit_qty(grid,P,R, p_fun, params)>0).*weight);
    
    %% Solve for price index
    price_puis = M * mu' * p_fun(grid).^(1-sig);
      
    %% compute the average revenu
    revenu = p_fun(grid).* profit_qty(grid,P,R, p_fun, params);
    rbar = mu' * revenu;
    
    
    %% Residual
    %RES(1)= R-L;
    RES(1)= P^(1-sig) - (price_puis);
    RES(2)= R -M*rbar;
     


end
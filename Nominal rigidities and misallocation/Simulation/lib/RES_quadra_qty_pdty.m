function RES=RES_quadra_qty_pdty(p, grid, params)

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
    

    phix=phi.*x;
    ptar= (w./(rho.*phi));
    %% eta'
    etaprim = xi./(phix).*sig.*(2/sig * (1./ptar-1./p)-(p./ptar-1).^2).*p.^(-1-sig); %Quadratique cost * quantity produce /pdty

    %% FOC
    RES = 1/rho * w ./ phix - etaprim/(sig-1).* p.^(1+sig)-p;
    

end
function ess= calc_ess(x)
    N = numel(x);
    
    [ac,lags] = xcov(x,'biased');
    ac = ac./ac(N);
    rho=ac(lags>0);
    
    max_pair=floor((numel(rho)-1)/2);
    rho_paired = zeros([max_pair 1]);
    
    for j = 1:max_pair
        rho_prop = rho(2*j)+rho(2*j+1);
        if rho_prop > 0
            rho_paired(j) = rho_prop;
        else
            rho_paired=rho_paired(1:j-1);
            break
        end
    end
    
    for j = 2:size(rho_paired,1)
        if rho_paired(j) > rho_paired(j-1)
            rho_paired(j) = rho_paired(j-1);
        end
    end
    
    P0 = 1+rho(1);
    ess = N / (-1+2*P0+2*sum(rho_paired));
end
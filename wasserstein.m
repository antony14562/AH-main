function ws = wasserstein(x,N,potential)
%define normalized pdf
    potential = lower(potential);
    
    switch potential
        case 'quartic'
            F = @(x) 0.25.*x.^4;
        case 'dwell_lip'
            F = @(x) ((x+1).^2 .* (x-1).^2) ./ (x.^2 + 1);
        case 'dwell_nonlip'
            F = @(x) 0.25.*x.^4-0.5.*x.^2;
        otherwise 
            error('unknown pot')
    end
       
    
    eF = @(x) exp(-2*F(x));
    Z = integral(eF,-Inf,Inf);
    norm_pdf = @(x) eF(x)./Z;
    
    %eval points
    xvals = linspace(-10,10,10^6);
    
    pdfvals = norm_pdf(xvals);  % evaluated pdf at each point
    cdfvals = cumtrapz(xvals,pdfvals); % cdf evaluation at each point
    
    [cdfvals,ia] = unique(cdfvals);
    xvals = xvals(ia);
    
    Finv = @(x) interp1(cdfvals, xvals, x); 
    
    %sample from interpolation
    U = unifrnd(0,1,([N 1]));
    samp = Finv(U);
    
    x = sort(x);
    q = sort(samp);
    ws = mean(abs(x-q));
    
   
    end


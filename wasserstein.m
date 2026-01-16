%define normalized pdf
F = @(x) 0.25.*x.^4-0.5.*x.^2;
%F = @(x) 0.25.*x.^4;
eF = @(x) exp(-2*F(x));
Z = integral(eF,-Inf,Inf);
norm_pdf = @(x) eF(x)./Z;

%eval points
xvals = linspace(-3,3,10^6);

pdfvals = norm_pdf(xvals);  % evaluated pdf at each point
cdfvals = cumtrapz(xvals,pdfvals); % cdf evaluation at each point

[cdfvals,ia] = unique(cdfvals);
xvals = xvals(ia);

Finv = @(x) interp1(cdfvals, xvals, x); 

%sample from interpolation
N = 10^5;
U = unifrnd(0,1,([N 1]));
samp = Finv(U);


function w = ws(x,q)
    x = sort(x);
    q = sort(q);
    w = mean(abs(x-q));
end

ws(samp,Q(:,2))


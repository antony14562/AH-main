%Use this for multiple runs with different deltas

N=8*3*10^6;
d=500;

sims=1;
beta = linspace(0.135,0.135,sims);%delta values
T=10;  
Delta = T/(d+1);

X = ones([N 3]);
acceptrates = zeros([sims 1]);
time_taken_avg = ones([sims 1]);

%F = @(x) 0.25.*x.^4; F1 = @(x) x.^3; F2 = @(x) 3.*x.^2; F3 = @(x) 6.*x; F4 = 6;

F = @(x) 0.25.*x.^4 - 0.5.*x.^2; F1 = @(x) x.^3 - x; F2 = @(x) 3.*x.^2 -1; F3 = @(x) 6.*x; F4 = 6;


%Construct L and cholesky of DeltaL = A for sampling
%L in this is -L from IMLA, see notes
l1 = -1*ones([d-1 1]); l2= 2*ones([d 1]);
L = (1/Delta^2) * spdiags([[l1;0] l2 [0;l1]], [-1 0 1], d, d);
R = chol(Delta.*L);

Phi = @(x) (Delta/2) * sum(F1(x).^2 - F2(x));

%% PCN


for s = 1:sims
    x0 = 0*ones([d 1]); %initial x
    accept_count=0;
    tic
    one_percent = tic;
    for j = 1:N
       
        xi = R \ normrnd(0,1,[d,1]); %solving 
        x1 = sqrt((1-beta(s)^2)).*x0 + beta(s).*xi;
        
        %accept/reject
        alpha = min(0,Phi(x0)-Phi(x1));
        r = unifrnd(0,1);
        if log(r) < alpha %accept
            x0 = x1;
            accept_count = accept_count +1;
        end
        
        %store first,mid, last components (points)
        
        X(j,1) = x0(1,:)';
        X(j,2) = x0(d/2,:)';
        X(j,3) = x0(end,:)';
        if mod(j,N/100) == 0  || j == N %tracker for if sim = 1
            one_percent_elapsed = toc(one_percent);
            fprintf('\r%d: %0.1f%%, elapsed time: %0.1f%',s,100*j/N,one_percent_elapsed) 
        end
    end
    time_taken_avg(s) = toc/N;
    acceptrates(s) = accept_count/N;
    disp(100*s/sims + "%")
end

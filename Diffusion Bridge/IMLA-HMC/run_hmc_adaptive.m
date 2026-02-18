function [Q,time,acceptrate,ws,ess,deltas] = run_hmc_adaptive(N,K,d,init_delta,target,T,pot,P,Uu)
%%%%%%%%%%%%
%N iterations
%K integration steps
%d dimension
%init_delta starting delta
%target target accept rate
%T bridge length
%pot potential "dwell_lip","dwell_nonlip"
%P random vars P randn(d, N)
%Uu random acceptance chance rand(1, N)
%%%%%%%%%%%%



% Select potential
if pot == "dwell_lip"
    U = @(Delta,L,x) (Delta/2)*x'*(-L)*x+Delta*sum(...
        -(x.^8 + 4*x.^6 + 18*x.^4 + 12*x.^2 - ...
        3 - 2*x.^10 -  8*x.^8 + 4*x.^6 + 24*x.^4 - 18*x.^2) ...
        ./((x.^2 + 1).^4));
    gradU = @(Delta,L,x) -L*x + ((4*x.*(x.^2 - 1).*(12*x.^2 + 12 + x.^8 ...
        + 6*x.^6 + 24*x.^4 + 42*x.^2 - 9))./((x.^2 + 1).^5));
    grad2U = @ (Delta,L,x) -L + spdiags( ...
        (4*(-60*x.^6 + 60*x.^4 + 108*x.^2 - 12 + ...
            x.^12 + 6*x.^10 - 9*x.^8 + 36*x.^6 + 447*x.^4 - 234*x.^2 + 9)) ...
        ./ ((x.^2 + 1).^6), ...
        0, numel(x), numel(x));
elseif pot == "dwell_nonlip"
    U = @(Delta,L,x) (Delta/2)*x'*(-L)*x+Delta*sum(x.^6/2 -x.^4 - x.^2 + 1/2);
    gradU = @(Delta,L,x) -L*x + (3*x.^5-4*x.^3-2*x);
    grad2U = @ (Delta,L,x) -L + spdiags( ( 15*x.^4 - 12*x.^2-2 ), ...
        0,numel(x),numel(x));
else
    error("Invalid potential. Select 'dwell_lip' or 'dwell_nonlip' as potential")
end

%%%%%%%%%%%%

tol = 10^-9; % Newton tol
Delta = T/(d+1);
coord = floor(d/2); % Midpoint

%Construct L
l1 = ones([d-1 1]); l2= -2*ones([d 1]);
L = (1/Delta^2) * spdiags([[l1;0] l2 [0;l1]], [-1 0 1], d, d);

%Construct Hess, Grad for newton
I = @(u,x,delta,Delta,z,L) u - x + delta * gradU(Delta,L,(0.5 .* (x + u))) ...
     - sqrt(2.*delta).*z;
Ip = @(u,x,delta,Delta,L) 0.5.*(delta)*grad2U(Delta,L,(0.5.*(x+u))) + speye(d);

%%%%%%%%%%%%%%%%%%%%%%%%% 
%% plotting vars
Q = zeros([N d]); %Final chain (storing whole chain)
%Q = zeros([N 3]) %Final chain (storing only 3 points)
iters = ones([N 1]); %Newton iteration counts
opt = ones([N 1]); %Newton final norm(update) val
deltas=NaN([N 1]); %delta values each iter
accepted = zeros([N 1]); %0 if iteration not accepted, 1 if accepted
accept_count = 0; %accepted count
%%%%%%%%%%%%%%%%%%%%%%%%% 

%%%%%%%%%%%%%%%%%%%%%%%%% 
%% Newton Function
function [q_n,iters,opt] = newton(q0,p,delta,tol,I,Delta,L,gradU,Ip)
    q_n = q0 - delta * gradU(Delta,L,q0)+ sqrt(2.*delta).*p;
    for jn = 1:10^4
        iters=jn;
        update = Ip(q_n,q0,delta,Delta,L)\I(q_n,q0,delta,Delta,p,L);
        q_n = q_n - update;
        if norm(update,Inf) < tol
            opt = norm(update,Inf);
            break
        end
        if jn == 10^4 disp("fail");opt = norm(update,Inf); end %max iterations
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%
%%  IMLA
init_q = 0*ones([d 1]);
q0 = init_q;%inital q
delta=init_delta;

tic;

for j = 1:N
    deltas(j) = delta;
    p0 = (1/sqrt(Delta))*P(:,j); 
    q=q0;p=p0;
    for k = 1:K %K HMC integration steps
        [q1,iters(j,1),opt(j,1)] = newton(q,p,delta,tol,I,Delta,L,gradU,Ip); %newton method to solve implicit step
        p1 = p - sqrt(2*delta).*gradU(Delta,L, (0.5*(q1 + q)) );
        q=q1;
        p=p1;
    end
    

    %accept/reject
    U1 = U(Delta,L,q); U0 = U(Delta,L,q0);
    K1 = 0.5*Delta*p'*p; K0 = 0.5*Delta*p0'*p0;
    alpha = min(0,U0-U1+K0-K1);
    if log(Uu(j)) < alpha % if accepted
        q0 = q;
        accept_count = accept_count +1;
        accepted(j) = 1;
    end
    delta=delta+delta*(accepted(j)-target)/j;
    %Save entire chain
    Q(j,:) = q0'; 
    
    %Save first,mid,last points
    %Must change ws calc line, Q construction line too.
    %Q(j,1) = q0(1,:)'; 
    %Q(j,2) = q0(coord,:)';
    %Q(j,3) = q0(end,:)';
end
Q = [init_q';Q]; %add initial q to start of chain var
time = toc;

acceptrate = accept_count/N;

%Change coord to 3 if saving only first,mid,last points
ws = wasserstein(Q(:,coord),N+1,pot);
ess = calc_ess(Q(:,coord));

end

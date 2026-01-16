%Use this for one run

N = 10^4;
d = 500;
tol = 10^-8;
delta = 0.0367; 
T = 10;
Delta = T/(d+1); %Delta = T/d+1

F = @(x) 0.25.*x.^4 - 0.5.*x.^2; F1 = @(x) x.^3 - x; F2 = @(x) 3.*x.^2 -1; F3 = @(x) 6.*x; F4 = 6;

%Construct L
l1 = ones([d-1 1]); l2= -2*ones([d 1]);
L = (1/Delta^2) * spdiags([[l1;0] l2 [0;l1]], [-1 0 1], d, d);

U = @(Delta,L,x) (Delta/2) * (-x'*L*x + sum(F1(x).^2 - F2(x)));
gradU = @(Delta,L,x) -Delta*L*x + Delta.*(F1(x).*F2(x) - 0.5*F3(x));
grad2U = @(Delta,L,x) -Delta*L + Delta*diag(F2(x).^2+F1(x).*F3(x)-0.5.*F4);

objF = @(u,x,delta,Delta,z) 2 * U(Delta,L,0.5 * (x + u)) + Delta./ (2 * delta) ...
            * norm(u - x - sqrt(2*delta/Delta).* z, 2);
I = @(u,x,delta,Delta,z,L) gradU(Delta,L,0.5 .* (x + u)) + ...
    Delta ./ delta .* (u - x - sqrt(2.*delta./Delta).* z);
Ip = @(u,x,delta,Delta,L) 0.5.*grad2U(Delta,L,0.5.*(x+u)) + (Delta/delta)*eye(d);

%plotting vars
Q = ones([N 3]);
iters = ones([N 1]);
opt = ones([N 1]);


%% Newton Function
function [q_n,iters,opt] = newton(q0,p,delta,tol,I,Delta,L,d,l2,F1,F2,F3,F4)
    q_n = q0;
    for j = 1:10^4
        iters=j;
        diag1= (Delta/delta) + 0.5.*((-1/Delta)*(l2)+Delta*(F2(0.5.*(q_n+q0)).^2+F1(0.5.*(q_n+q0)).*F3(0.5.*(q_n+q0))-0.5.*F4));
        diag2 = 0.5.*((-1/Delta)*ones([d-1 1]));
        q_n = q_n - tridiagonal_vector(diag2,diag1,diag2,I(q_n,q0,delta,Delta,p,L));
        if norm(I(q_n,q0,delta,Delta,p,L),Inf) < tol
            opt = norm(I(q_n,q0,delta,Delta,p,L),Inf);
            break
        end
        if j == 10^4-1 disp("fail"); end %max iterations
    end
end

%% IMLA 

tic;
one_percent = tic; % tracker for elapsed time track

accept_count = 0;
q0 = 0*ones([d 1]); %inital q

for j = 1:N
    p = normrnd(0,1,([d 1]));
    [q_h,iter_count,result] = newton(q0,p,delta,tol,I,Delta,L,d,l2,F1,F2,F3,F4); %newton method to solve implicit step
    iters(j,1) = iter_count;
    opt(j,1) = result;
    q1 = q_h;
    p1 = p - sqrt(2*delta/Delta).*gradU(Delta,L,q1/2 + q0/2); %p_n
    p1=-p1;

    %accept/reject
    U1 = U(Delta,L,q1); U0 = U(Delta,L,q0);
    K1 = 0.5*p1'*p1; K0 = 0.5*p'*p;
    alpha = min(0,-U1 + U0 - K1 + K0);
    r = unifrnd(0,1);
    if log(r) < alpha %accept
        q0 = q1;
        accept_count = accept_count +1;
    end

    %store only first,mid,last components(points)

    Q(j,1) = q0(1,:)';
    Q(j,2) = q0(d/2,:)';
    Q(j,3) = q0(end,:)';

    if mod(j,N/100) == 0  || j == N %tracker for if sim = 1
        one_percent_elapsed = toc(one_percent);
        fprintf('\n%0.1f%%, elapsed time: %0.1f%',100*j/N,one_percent_elapsed)
    end

end

time = toc;
acceptrates = accept_count/N;

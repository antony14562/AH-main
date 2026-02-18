clear all

%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parameters 
T_vals = [10];
delta_vals = [0.001,0.005:0.005:0.1];
d_vals = [20,100];
K_vals = [1];
Delta_vals = T_vals(:) ./ (d_vals(:).' + 1);
N=10^5;

nTests = numel(T_vals)*numel(delta_vals)*numel(d_vals);


c=0;
for j_K = 1:numel(K_vals)
    K = K_vals(j_K);
    for j_T = 1:numel(T_vals)
        T = T_vals(j_T);
        for j_del = 1:numel(delta_vals)
            delta = delta_vals(j_del);
            for j_dim = 1:numel(d_vals)
                d = d_vals(j_dim);
                rng(123)
                P = randn(d, N);
                Uu = rand(1, N);
                [Q,time,acceptrate,ws,ess] = run_hmc(N,K,d,delta,T,pot,P,Uu);
                
            end
        end
    end
end
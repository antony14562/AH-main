clear all

%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parameters 
T_vals = [10,50];
delta_vals = [0.03,0.08,0.12];
d_vals = [10:10:100,120:40:500,600,750,1000,2000];
K_vals = [1];

% T_vals = [10,20];
% delta_vals = [0.01,0.05];
% d_vals = [10,20,30];
% K_vals = [1];

Delta_vals = T_vals(:) ./ (d_vals(:).' + 1);
N=5*10^4;

nTests = numel(T_vals)*numel(delta_vals)*numel(d_vals);


c=0;
for j_K = 1:numel(K_vals)
    K = K_vals(j_K);
    for j_del = 1:numel(delta_vals)
        delta = delta_vals(j_del);
        for j_T = 1:numel(T_vals)
            T = T_vals(j_T);
            acceptrates = NaN([numel(d_vals) 1]);
            for j_dim = 1:numel(d_vals)
                d = d_vals(j_dim);
                c=c+1;
                rng(123)
                P = randn(d, N);
                Uu = rand(1, N);

                fprintf('\rd=%d,T=%d,delta=%0.3f,K=%d,test %d/%d',d,T,delta,K,c,nTests)
                
                [Q,time,acceptrates(j_dim),ws,ess] = run_hmc(N,K,d,delta,T,'dwell_nonlip',P,Uu);                
            end
            figure(j_del)
            plot(log(Delta_vals(j_T,:)),acceptrates,'LineStyle','-','Marker','o', ...
                'DisplayName', sprintf('T = %g', T));
            title('\delta = delta');xlabel('dim');ylabel('accept rate')
            legend(gca,'show','Location','best');
            hold on;
        end
    end
end
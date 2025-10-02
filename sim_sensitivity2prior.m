function sim_sensitivity2prior

run_analysis(10, 1);
run_analysis(10, 1000);
run_analysis(20, 1000);

end

function run_analysis(N, c)

K = 20;
T = 100;
num_sim = 1e3;
seed = 0;
pb = 0.05;

% -------------------------------------------------------------------------
[ep_thr] = run_null(N, K, .05, c);

rng(seed);
r = zeros(N, num_sim, K);
loglik = zeros(N, num_sim, K);            

mtrue = repmat([N zeros(1, K-1)], num_sim, 1);
for k=1:K
    p = zeros(N, T, num_sim);
    for i=1:num_sim    
        ph = pb + (1-pb)*rand(mtrue(i,k), T);
        pc = rand(N-mtrue(i,k), T);
        p(:, :, i) = cat(1, ph, pc);        
    end
    logp = squeeze(sum(log(p+eps), 2));    
    loglik(:, :, k) = logp;                
    
    r(:, :, k) = exp(logp);      
end

m = r./sum(r, 3);
m = log(squeeze(sum(m, 1)));
posterior = exp(m)+c;

[ep] = compute_exceedance(posterior, 1e4);

power = sum(ep(:, 1)>ep_thr);

fprintf('N=%d, alpha=%d, power=%d/1000\n', N, c, power);
end


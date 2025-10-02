function winner = ffx_null()
N = 50;
K = 3;
T = 200;
num_sim = 1e4;
log_bf_crit = 3;
seed = 0;
prior_dirichlet = 1;
num_samples_4ep = 1e4;

% -------------------------------------------------------------------------
num_batch = 1;
fname = fullfile('sum', sprintf('%s_logBF%d.mat', mfilename, log_bf_crit));
winner = struct('fixed', [], 'random', []);
if ~exist(fname, 'file')

    [ep_thr] = run_null(N, K);
    
    num_sim = num_sim/num_batch;
    rng(seed);
    for batch = 1:num_batch
        f = zeros(num_sim, K);
        r = zeros(N, num_sim, K);
        for k=1:K            
            p = rand(N, T, num_sim);

            logp = squeeze(sum(log(p+eps), 2));
            f(:, k) = sum(logp, 1);

            logp = squeeze(sum(log(p+eps), 2));            
            
            r(:, :, k) = exp(logp);
        end
        f = f - max(f, [], 2);

        r = r./sum(r, 3);
        r = (squeeze(sum(r, 1)));

        % posterior Dirichlet parameter
        posterior = r + prior_dirichlet; 
        
        % calculate exceedance probability
        [ep] = compute_exceedance(posterior, num_samples_4ep);

        for k=1:K
            kk = 1:K; kk(k) = [];
            f_evidence = f(:, k) - f(:, kk);
            f_evidence = min(f_evidence, [], 2);
            winner_batch.fixed(batch, k) = mean(f_evidence > log_bf_crit);

            winner_batch.random(batch, k) = mean(ep(:, k) > ep_thr);
        end
    end
    winner.fixed = mean(winner_batch.fixed, 1);
    winner.random = mean(winner_batch.random, 1);
    
    config = table2struct(table(N, K, T, num_sim, log_bf_crit, seed));    
    save(fname, 'config', 'winner_batch', 'winner');    
end
    
f = load(fname);
winner = f.winner;

end




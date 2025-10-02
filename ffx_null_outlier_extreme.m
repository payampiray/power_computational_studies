function winner = ffx_null_outlier_extreme()
N = 50;
K = 3;
T = 200;
p_min = 0;
p_max = 1;
n_outlier = 1;
p_outlier = .9;
num_sim = 1e4;
log_bf_crit = 3;
seed = 0;

% -------------------------------------------------------------------------
num_batch = 1;
fname = fullfile('sum', sprintf('%s_logBF%d.mat', mfilename, log_bf_crit));
winner = struct('fixed', [], 'random', []);
if ~exist(fname, 'file')
    
    num_sim = num_sim/num_batch;
    rng(seed);
    for batch = 1:num_batch
        f = zeros(num_sim, K);
        r = zeros(N, num_sim, K);

        for k=1:K            
            p = p_min + (p_max - p_min)*rand((N-n_outlier), T, num_sim);
            if k>1
                po = p_min + (p_max - p_min)*rand(n_outlier,T,num_sim);
            else
                po = p_outlier + (1-p_outlier)*rand(n_outlier,T, num_sim);
            end
            p = cat(1, p, po);

            logp = squeeze(sum(log(p+eps), 2));
            f(:, k) = sum(logp, 1);

        end
        f = f - max(f, [], 2);

        for k=1:K
            kk = 1:K; kk(k) = [];
            f_evidence = f(:, k) - f(:, kk);
            f_evidence = min(f_evidence, [], 2);
            winner_batch.fixed(batch, k) = mean(f_evidence > log_bf_crit);

        end
    end
    winner.fixed = mean(winner_batch.fixed, 1);
    
    config = table2struct(table(N, K, T, n_outlier, p_outlier, num_sim, log_bf_crit, seed));    
    save(fname, 'config', 'winner');    
end

f = load(fname);
winner = f.winner;    




end




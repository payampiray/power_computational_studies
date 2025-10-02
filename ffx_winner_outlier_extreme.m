function clear_win = ffx_winner_outlier_extreme()
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
if ~exist(fname, 'file')
    
    num_sim = num_sim/num_batch;
    rng(seed);

    for batch = 1:num_batch    
        f2 = zeros(0, K);
        r2 = zeros(0, K);
        num_win2 = 0;
        while num_win2 < num_sim
            f = zeros(num_sim, K);
            r = zeros(N-n_outlier, num_sim, K);
            for k=1:K
                p = p_min + (p_max - p_min)*rand((N-n_outlier), T, num_sim);
                logp = squeeze(sum(log(p+eps), 2));
                f(:, k) = sum(logp, 1);
                r(:, :, k) = exp(logp);
            end
            f2_win = f(:, 2) - f;
            win2 = sum(f2_win > log_bf_crit, 2) == 2;
            
            f2 = cat(1, f2, f(win2, :));
            num_win2 = size(f2, 1);

            r = r./sum(r, 3);
            r = (squeeze(sum(r, 1)));   

            r2 = cat(1, r2, r(win2, :));
        end
        f = f2(1:num_sim, :);
        
        num_sim = size(f, 1);
        for k=1:3
            if k>1
                po = p_min + (p_max - p_min)*rand(n_outlier, T, num_sim);
            else
                po = p_outlier + (1-p_outlier)*rand(n_outlier, T, num_sim);
            end
            logp = zeros(n_outlier, num_sim);
            logp(1:n_outlier, :) = (sum(log(po+eps), 2));
        
            fo = sum(logp, 1)';
            f(:, k) = f(:, k) + fo;

        end

            
        for k=1:K
            kk = 1:K; kk(k) = [];
            f_evidence = f(:, k) - f(:, kk);
            f_evidence = min(f_evidence, [], 2);
            clear_win_batch.fixed(batch, k) = mean(f_evidence > log_bf_crit);
        end
    end
    clear_win.fixed = mean(clear_win_batch.fixed, 1);

    config = table2struct(table(N, K, T, n_outlier, p_outlier, num_sim, log_bf_crit, seed));
    save(fname, 'config', 'clear_win_batch', 'clear_win');    
end


f = load(fname);
clear_win = f.clear_win;

end

function [N, K, T, p_min, p_max, n_outlier, p_outlier, num_sim, log_bf_crit, seed] = untitled
N = 50;
K = 3;
T = 200;
p_min = 0;
p_max = 1;
n_outlier = 1;
p_outlier = .9;
num_sim = 10^4;
log_bf_crit = 3;
seed = 0;
end
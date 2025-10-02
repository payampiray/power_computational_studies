function winner = ffx_null_outlier_modest()
N = 50;
K = 3;
T = 200;
p_min = 0;
p_max = 1;
p_outlier = .1;
num_sim = 1e4;
log_bf_crit = 3;
seed = 0;

% -------------------------------------------------------------------------
ffx = [];
for n_outliers = 0:5
    n_outlier = n_outliers;
    fname = fullfile('sum', sprintf('%s_%doutliers.mat', mfilename, n_outlier));
    if ~exist(fname, 'file')
        
        rng(seed);
        f = zeros(num_sim, K);    
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

        winner = zeros(1, K);
        for k=1:K
            kk = 1:K; kk(k) = [];
            f_evidence = f(:, k) - f(:, kk);
            f_evidence = min(f_evidence, [], 2);
            winner(k) = mean(f_evidence > log_bf_crit);
        end

        config = table2struct(table(N, K, T, n_outlier, p_outlier, num_sim, log_bf_crit, seed));    
        save(fname, 'config', 'winner');    
    end
    
    loaded = load(fname);
    ffx = cat(1, ffx, loaded.winner);    
end

%--------------------------------------------------------------------------
% close all;

x = 0:(size(ffx, 1)-1);
y = ffx*100;
`
% T = array2table([x' y], 'VariableNames', {'Number of outliers', 'Model 1', 'Model 2', 'Model 3'});
% writetable(T,'source_Fig3.csv');

labels = {'Model 1', 'Model 2', 'Model 3'};


fs = 14;
fsy = 16;

fsiz = [0 0 .3 .3];
figure; set(gcf,'units','normalized'); set(gcf,'position',fsiz);

bar(x, y);
set(gca, 'ylim', [0 100], 'FontSize', fs, 'xtick', x);

legend(labels, 'FontSize', fsy, 'Location','northwest', 'orientation', 'vertical', 'box', 'off', 'AutoUpdate', 'off');

ylabel('Winner of model selection', 'fontsize', fsy);
xlabel('Number of outliers (out of 50)', 'fontsize', fsy);

set(gca, 'box', 'off', 'ticklength', [0 0 ]);

end




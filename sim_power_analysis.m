function sim_power_analysis

sample_sizes = 20:30:200;
model_space = [4 7 10];

clc;
fname = fullfile('sum', sprintf('%s.mat', mfilename));
if ~exist(fname, 'file')
    for n=1:length(sample_sizes)
        N = sample_sizes(n);
        power = nan(1, length(model_space));
        for k=1:length(model_space)
            K = model_space(k); 
            [exceedance_threshold] = run_null(N, K);
            [power(k)] = run_power(N, K, exceedance_threshold);
            labels{k} = sprintf('K%d', K);
            fprintf('N=%02d, K=%02d is done| power=%0.2f\n', N, K, power(k));
        end
        powers(n, :) = power;        
    end

    T = array2table([sample_sizes' powers], 'VariableNames', ['sample size', labels]);
    save(fname, 'T', 'model_space');
end
f = load(fname);
T = f.T;
model_space = f.model_space;

labels = cellstr(num2str(model_space'));
powers = table2array(T(:, 2:4));
powers = round(powers*100)/100;


% T = array2table(table2array(T), 'VariableNames', {'Sample size', 'Model space size of 4', 'Model space size of 7', 'Model space size of 10'});
% writetable(T,'source_Fig1.csv');

%--------------------------------------------------------------------------

x = sample_sizes;
y = powers;


fs = 14;
fsy = 18;
col = [0    0.4470    0.7410; 0.8500    0.3250    0.0980];

fsiz = [0 0 .3 .3];
h = figure; set(gcf,'units','normalized'); set(gcf,'position',fsiz);

bar(x, y);
set(gca, 'ylim', [0 1], 'FontSize', fs, 'xtick', x);
% % plot(x, y, 'linewidth', 2);

hg = legend(labels, 'FontSize', fsy, 'Location','northwest', 'orientation', 'horizontal', 'box', 'off', 'AutoUpdate', 'off');
title(hg, 'Model space size', 'FontWeight','normal');

ylabel('Power', 'fontsize', fsy);
xlabel('Sample size', 'fontsize', fsy);

xl = get(gca, 'xlim');

% hold on;
% plot(xl, .8*[1 1], 'linewidth', 2, 'Color', 'k');

set(gca, 'box', 'off', 'ygrid', 'on', 'ticklength', [0 0 ]);

end
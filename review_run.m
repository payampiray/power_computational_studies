function review_run()

fname = fullfile('sum', sprintf('%s.mat', mfilename));

if ~exist(fname, 'file')
    review = load(fullfile('sum', 'review_run.mat')); review = review.review;

    studies = review(:, {'sample_size', 'model_size', 'sample_size_extra', 'model_size_extra'});
    idx_to_compute_power = 1:size(studies.sample_size, 1);
    [studies] = compute_study_power(studies, idx_to_compute_power);

    review.thr = studies.thr;
    review.false_positives = studies.false_positives;
    review.power = studies.power;
    review.thr_extra = studies.thr_extra;
    review.false_positives_extra = studies.false_positives_extra;
    review.power_extra = studies.power_extra;
    save(fname, 'review');    
end
f = load(fname); review = f.review;

power1 = review.power;
power_extra = review.power_extra;
power = nanmean([power1, power_extra], 2);

power = ceil(power*100)/100;

required_power = .8;
low_power = sum(power<required_power);

[~, idx] = sort(power, 'descend');
[sorted_power] = power(idx);

sorted_fixed = review.method(idx);
sorted_fixed = strcmp(sorted_fixed, 'fixed');

% writematrix(sorted_power,'source_Fig2.csv');

%--------------------------------------------------------------------------
fs = 12;
fsy = 16;

colmap = [201 92 46; 228 179 69]/255;

fsiz = [0 0 .5 1];
figure; set(gcf,'units','normalized'); set(gcf,'position',fsiz);


K = length(sorted_power);
wb = .6;
for i=1:K
    col = colmap(1, :);
    if sorted_fixed(i), col = colmap(2, :); end
    
    hb(i) = barh(i, sorted_power(i),wb,'FaceColor',col); hold on;

end

set(gca, 'fontsize', fs);
ylabel('Publication number', 'FontSize', fsy);
xlabel('Estimated power', 'FontSize', fsy);

yl = get(gca, 'ylim');

ytick = 1:3:K;

hold on;
plot(required_power*[1 1], yl, 'linewidth', 2, 'color', 'k');

set(gca, 'box', 'off', 'ytick', ytick, 'xgrid', 'on');
end


function [studies] = compute_study_power(studies, idx_to_compute_power, seed)
if nargin<3, seed = 0; end

sample_size = studies.sample_size;
model_size = studies.model_size;
sample_size_extra = studies.sample_size_extra;
model_size_extra = studies.model_size_extra;

size_extra = nan(size(sample_size, 1), 1);
for i= 1:size(sample_size, 1)
    size_extra(i) = size(sample_size_extra{i}, 2);
end

power = nan(size(model_size));
thr = nan(size(model_size));
false_positives = nan(size(model_size));

thr_extra = nan(size(model_size, 1), max(size_extra));
false_positives_extra = nan(size(model_size, 1), max(size_extra));
power_extra = nan(size(model_size, 1), max(size_extra));

for l = 1:length(idx_to_compute_power)        
    i = idx_to_compute_power(l);
    
    [thr(i), false_positives(i)] = run_null(sample_size(i), model_size(i), seed);
    [power(i)] = run_power(sample_size(i), model_size(i), thr(i), seed);
    
    if ~isempty(sample_size_extra{i})
        for j=1:size(sample_size_extra{i}, 2)
            [thr_extra(i, j), false_positives_extra(i, j)] = run_null(sample_size_extra{i}(j), model_size_extra{i}(j), seed);
            [power_extra(i, j)] = run_power(sample_size_extra{i}(j), model_size_extra{i}(j), thr_extra(i, j), seed);
        end
    end
    
    fprintf('%02d/%02d is done\n', l, size(idx_to_compute_power, 1));
end

studies.thr = thr;
studies.false_positives = false_positives;
studies.power = power;
studies.thr_extra = thr_extra;
studies.false_positives_extra = false_positives_extra;
studies.power_extra = power_extra;

end
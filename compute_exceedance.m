function [exceedance] = compute_exceedance(dirichlet_parameters, num_samples_4ep)
if nargin<2
    num_samples_4ep = 1e4;
end

num_samples_4ep_batch = 1e3;
num_batch = ceil(num_samples_4ep/num_samples_4ep_batch);
sum_max = [];
for i=1:num_batch
    posterior_samples_batch = dirichletrnd3x(dirichlet_parameters, num_samples_4ep_batch);
    [~, sum_max_batch] = max(posterior_samples_batch, [], 2);
    sum_max_batch = squeeze(sum_max_batch);
    sum_max = cat(2, sum_max, sum_max_batch);
end

exceedance = nan(size(dirichlet_parameters));
for k=1:size(exceedance, 2)
    exceedance(:, k) = mean(sum_max == k, 2);
end

end



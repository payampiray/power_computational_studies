function [power] = run_power(num_participants, num_models, exceedance_threshold, seed)
if nargin<4
    seed = 0;
end

if ~exist('temp', 'dir'), mkdir('temp'); end
fname = fullfile('temp',sprintf('%s_N%dK%d.mat', 'power', num_participants, num_models));

if ~exist(fname, 'file')    
    num_sim = 1e4;
    prior_parameter = 1;
    max_iter = 100;
    num_samples_4ep = 1e4;
    
    rng(seed);
    
    % generate true population-level scenarios
    L = 0;
    population_samples = [];
    num_iter = 0;
    while (L<num_sim)&&(num_iter < max_iter)
        r = dirichletrnd(num_sim, ones(1,num_models));
        model1_best_idx = (r(:, 1)>max(r(:, 2:end), [], 2));
        population_samples = cat(1, population_samples, r(model1_best_idx, :));
    
        L = length(population_samples);
        num_iter = num_iter+1;
    end
    population_samples = population_samples(1:num_sim, :);
    
    % generate group samples (multinomial distribution)
    group_samples = mnrnd(num_participants, population_samples);
    
    % posterior Dirichlet parameters
    posterior_parameters = group_samples + prior_parameter; 
    
    % calculate exceedance probabilities based on posterior Dirichlet parameters
    [exceedance_prob] = compute_exceedance(posterior_parameters, num_samples_4ep);
    
    config = struct('num_sim', num_sim, 'max_iter', max_iter, 'num_samples_4ep', num_samples_4ep, 'prior_parameter', prior_parameter, 'seed', seed);
    save(fname, 'exceedance_prob', 'posterior_parameters', 'population_samples', 'config');
end

f = load(fname);
exceedance_prob1 = f.exceedance_prob(:, 1);
power = mean(exceedance_prob1 > exceedance_threshold);


end


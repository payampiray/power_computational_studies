function [power]=run_effectsize(num_participants, num_models, exceedance_threshold, seed)
if nargin<4
    seed = 0;
end

if ~exist('temp', 'dir'), mkdir('temp'); end
fname = fullfile('temp',sprintf('%s_K%d_seed%d.mat', 'effectsize', num_models, seed));
if ~exist(fname, 'file')      
    tol_effectsize = .01;
    step_base = 0.02;
    true_dirichlet_parameter = 1;
    effectsize_acceptable = 0.3;

    [effectsize] = compute_effect_size(num_models, true_dirichlet_parameter, seed);
    effectsizes(1) = mean(effectsize);
    
    true_dirichlet_parameters(1) = true_dirichlet_parameter;

    if effectsizes< effectsize_acceptable
        step_size = -step_base;
    else
        step_size = +step_base;
    end

    i = 1;
    max_iter = 50;
    def_pre = effectsizes(1) - effectsize_acceptable;
    while (abs(effectsizes(end) - effectsize_acceptable)>tol_effectsize)&&(i < max_iter)
        def = effectsizes(end) - effectsize_acceptable;
        if (def>0) ~= (def_pre>0)
            step_size = -step_size/2;
        end

        i = i+1;
        parameter_before = true_dirichlet_parameter;
        true_dirichlet_parameter = parameter_before + step_size;
        true_dirichlet_parameters(i) = true_dirichlet_parameter;
        [effectsize] = compute_effect_size(num_models, true_dirichlet_parameter, seed);
        effectsizes(i) = mean(effectsize);

        def_pre = def;
    end
    

    config = struct('effectsize_acceptable', effectsize_acceptable, ...
                    'step_base', step_base, 'tol_effectsize', tol_effectsize, 'seed', seed);
        
    effect = struct('true_dirichlet_parameter', true_dirichlet_parameter, 'effect_size', effectsizes(end));
    save(fname, 'effect', 'config');    
end

f = load(fname);
effect = f.effect;
true_dirichlet_parameter = effect.true_dirichlet_parameter;


fname = fullfile('temp',sprintf('%s_N%dK%d_seed%d.mat', 'power_effectsize', num_participants, num_models, seed));
if ~exist(fname, 'file')     
    % run power
    prior_parameter = 1;
    num_samples_4ep = 1e4;

    rng(seed);
        
    [~, population_samples] = compute_effect_size(num_models, true_dirichlet_parameter, seed);
    % generate group samples (multinomial distribution)
    group_samples = mnrnd(num_participants, population_samples);
    
    % posterior Dirichlet parameters
    posterior_parameters = group_samples + prior_parameter;
    
    % calculate exceedance probabilities based on posterior Dirichlet parameters
    [exceedance_prob] = compute_exceedance(posterior_parameters, num_samples_4ep);
        
    config = struct('num_samples_4ep', num_samples_4ep, 'prior_parameter', prior_parameter, 'seed', seed);
    save(fname, 'exceedance_prob', 'posterior_parameters', 'population_samples', 'effect', 'config');
end


f = load(fname);
exceedance_prob = f.exceedance_prob;
exceedance_prob1 = exceedance_prob(:, 1);

power = mean(exceedance_prob1 > exceedance_threshold);
end

function [effectsize, population_samples] = compute_effect_size(num_models, prior, seed)
max_iter = 500;
num_sim  = 1e4;

rng(seed);

% generate true population-level scenarios
L = 0;
population_samples = []; effectsize = [];
num_iter = 0;
while (L<num_sim)&&(num_iter < max_iter)
    r = dirichletrnd(num_sim, prior*ones(1,num_models));
    model1_best_idx = (r(:, 1)>max(r(:, 2:end), [], 2));
    diff1_2nd = (r(:, 1) - max(r(:, 2:end), [], 2));
    effectsize = cat(1, effectsize, diff1_2nd(model1_best_idx, :));
    population_samples = cat(1, population_samples, r(model1_best_idx, :));

    L = length(population_samples);
    num_iter = num_iter+1;
end
effectsize = effectsize(1: num_sim, :);
population_samples = population_samples(1: num_sim, :);

end


function [exceedance_threshold, false_positives_rate] = run_null(num_participants, num_models, seed, alpha_acceptable, prior_parameter, num_sim, num_samples_4ep)

if nargin<3
    seed = 0;
end
if nargin<4
    alpha_acceptable = 0.05;
end

if ~exist('temp', 'dir'), mkdir('temp'); end
fname = fullfile('temp',sprintf('%s_N%dK%d_seed%d.mat', 'null', num_participants, num_models, seed));
if ~exist(fname, 'file')
    
    if nargin<5
        prior_parameter = 1;
    end
    if nargin<6
        num_sim = 1e4;
    end
    if nargin<7
        num_samples_4ep = 1e4;
    end    
        
    rng(seed);
    
    group_samples = mnrnd(num_participants, ones(num_sim, num_models)/num_models);
    posterior_parameters = group_samples + prior_parameter;
    
    [exceedance_prob] = compute_exceedance(posterior_parameters, num_samples_4ep);
    
    [exceedance_threshold, false_positives_rate] = find_critical_threshold(exceedance_prob, alpha_acceptable);  

    config = struct('num_sim', num_sim, 'num_samples_4ep', num_samples_4ep, 'prior_parameter', prior_parameter, ...
                    'alpha_acceptable', alpha_acceptable, 'seed', seed);
    save(fname, 'posterior_parameters', 'exceedance_prob', 'exceedance_threshold', 'false_positives_rate', 'config');
end

f = load(fname);
exceedance_threshold = f.exceedance_threshold;
false_positives_rate = [];
if nargout>1
    false_positives_rate = f.false_positives_rate;
end

end


function [exceedance_threshold, false_positives_rate] = find_critical_threshold(exceedance_prob, alpha_acceptable)
threshold_min = 1/size(exceedance_prob, 2);
threshold_max = 1; 
step_level = 1e-4;

num_samples = size(exceedance_prob, 1);
exceedance_prob_max = max(exceedance_prob, [], 2);

thresholds = threshold_min: step_level: threshold_max;

false_positives_rate = sum(exceedance_prob_max>=thresholds)/num_samples;

idx = find( false_positives_rate <= alpha_acceptable, 1, 'first');
exceedance_threshold = thresholds(idx);
false_positives_rate = false_positives_rate(idx);

end

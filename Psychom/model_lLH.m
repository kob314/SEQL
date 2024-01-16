function lPost = model_lLH(par_v,par_info,D)

%% transform parameters
par = parameter_wrap(par_v,par_info,'v2s');

target = 0;

%% priors
% global
    % sigma
    target  = target + normpdf(par.bias_subj_sigma_tf,2,1);
    par.bias_subj_sigma = exp(-par.bias_subj_sigma_tf);
    % kappa
    target  = target + sum(log(normpdf(par.kappa,0,2)));
    % tau
    target  = target + sum(log(normpdf(par.tau_tf,1,1)));
    par.tau = exp(-par.tau_tf);

% target = target + sum(log(normpdf(par.bias_exp,0,2)));
% target = target + sum(log(exppdf(par.bias_exp,0,2)));

% experiment-wise
    % bia                                                                                                                                                                                                                                                                                                                          s_exp
    target = target + sum(log(normpdf(par.bias_exp,0,2)));

% subject-wise
    % bias_subj
    target = target + sum(log(normpdf(par.bias_subj,0,par.bias_subj_sigma)));
    % stim
    target = target + sum(log(normpdf(par.stim,0,2)));
    % lambda 
    target  = target + sum(log(normpdf(par.lambda_tf,-2,1)));
    par.lambda = 1./(1+exp(-par.lambda_tf));


%% likelihood
w_ST = exp(-par.tau*[1:size(D.past_r,1)]);
w_ST = w_ST./sum(w_ST,2);
STSE = w_ST * D.past_r;
% size(w_ST)
% size(STSE)


mu    = par.bias_exp(D.exp) + par.bias_subj(D.subj) + par.stim(D.subj).*(D.z-.5)*2.*D.y + par.kappa*STSE;
mu_tf = 1 ./ ( 1 + exp(-mu) );

P_r1 = mu_tf.*(1-par.lambda(D.subj))+.5*par.lambda(D.subj);

lLH = sum(D.r .* log(P_r1) + (1-D.r) .* log(1-P_r1),2);

%% posterior
lPost = lLH + target; 


% mu = beta_bias_subjFH[subj_idxs]*FH + beta_bias_expFH[exp_idxs]*FH + beta_bias_subjSH[subj_idxs]*(1-FH) + beta_bias_expSH[exp_idxs]*(1-FH) + beta_stim[subj_idxs]*(z-.5)*2 * y + beta_past_kappa*STSE 

% P_r1_G_tr = mu_tf*(1-par.lambda(1))+.5*par.lambda(1);
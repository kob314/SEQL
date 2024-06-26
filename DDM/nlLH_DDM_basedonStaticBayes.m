function [nlLH p] = nlLH_DDM_basedonStaticBayes(par_v,par_info,data_in,o)

% par_v
par = parameter_wrap(par_v,par_info,'v2s');

if ~isfield(par,'s')
    w_bias =  .5;
    v_bias =   0;
elseif isfield(par,'wb')
    % %% bias transformation
    % bn = data_in.bias_M + par.s;
    bn = sign(data_in.bias_M).*abs(data_in.bias_M) .^ par.s;
    % %% v and w computation
    % w_bias =  par.wb * bn +.5;
    % v_bias =  par.vb * bn;
    w_bias =  par.wb * bn +.5;
    v_bias =  sign(par.wb)*par.vb * bn;
    
else
    %% AP(q)-RV(nu) transformation
    qn  = data_in.q_M   + data_in.shift*par.s;
    nun = data_in.nu_M  - data_in.shift*par.s;
    %% v and w computation
    w_bias =  par.w * (qn-.5) + .5;
    v_bias = -par.v * (nun-.5);
end

if ~isfield(par,'kappa_w')
    par.kappa_w = 0;
    par.kappa_v = 0;
end

w_aux = w_bias + par.kappa_w*data_in.STSE;
w     = 1./(1+exp(-4*(w_aux-.5)));
v_aux = v_bias + par.v_amp*data_in.yz + par.kappa_v*data_in.STSE;
v     = par.a*v_aux;


%% probability distributions
p = wfpt_t0noise_vec(reshape(v,1,1,[]),reshape(w,1,1,[]),par.a,par.t0,par.t0_lsig,par.l,o.grid_t,o.k);


nlLH = -sum(log(p(data_in.p_idx)));

if isinf(nlLH)
    nlLH = 10^10;
end




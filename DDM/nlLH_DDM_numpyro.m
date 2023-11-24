function [nlLH p] = nlLH_DDM_numpyro(par_v,par_info,data_in,o)

% par = parameter_wrap(par_v,par_info,'v2s');
par = par_v;

par.t0_lsig     = (par.t0_lsig_aux+.01)/par.t0;

stim_strength = 2*(data_in.z-.5).*data_in.y;



qn  = par.q_exp(data_in.e);
nun = o.true_AP(data_in.e) - par.q_exp(data_in.e) +.5;

% rrt

% scatter(data_in.q_M,data_in.nu_M)
% hold on
% scatter(qn,nun)
% plot([0 1; .5 .5]',[.5 .5;0 1]','k')
% plot([-1 1; -1 1]',[[1 -1]+1;[1 -1]+1.15]','k')
% rrt


w_bias = par.ws * (qn-.5) + .5;  %par.w *  (sign( data_in.q_M-.5  + par.s).*abs( data_in.q_M-.5  + par.s).^par.sw) +.5;  % q dependent
v_bias = -par.vs * (nun-.5);%par.v * -(sign(data_in.nu_M-.5  - par.s).*abs(data_in.nu_M-.5  - par.s).^par.sv);      % nu dependent


% subplot(1,2,1)
% scatter(data_in.q_M,w_bias)
% subplot(1,2,2)
% scatter(data_in.nu_M,v_bias)


% w_bias = par.w * (data_in.q_M(data_in.e)-.5    + par.s) +.5;  % q dependent
% v_bias = par.v * -(data_in.nu_M(data_in.e)-.5  - par.s);    % nu dependent

% v_bias = par.v * -(data_in.nu_M(data_in.e)-.5 );
% w_bias = par.w * (data_in.q_M(data_in.e)-.5 ) +.5;

if ~isfield(par,'kappa_w')
    par.kappa_w = 0;
    par.kappa_v = 0;
    par.tau_w   = 0;
end

weight_ST  = exp(-par.tau_w*[1:o.n_back]);
weight_ST  = weight_ST/sum(weight_ST);
past_r = (data_in.past_r-.5)*2;
past_r(isnan(past_r)) = 0;
data_in.STSE = weight_ST*past_r;


w     = 1./(1+exp(-4*(w_bias-.5 + par.kappa_w*data_in.STSE)));
v     = par.a*(v_bias + par.v_amp*stim_strength + par.kappa_v*data_in.STSE);




p = wfpt_t0noise_vec(reshape(v,1,1,[]),reshape(w,1,1,[]),par.a,par.t0,par.t0_lsig,par.l,o.grid_t,o.k);

nlLH = -sum(log(p(data_in.p_idx)));



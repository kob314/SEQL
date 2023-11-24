function [nlLH p] = nlLH_DDM_basedonStaticBayes_restricted(par_v,par_info,data_in,o)

par = parameter_wrap(par_v,par_info,'v2s');

if isfield(par,"s")
    bn = data_in.bias_M + par.s;
    w_bias =  par.wb * bn +.5;  %par.w *  (sign( data_in.q_M-.5  + par.s).*abs( data_in.q_M-.5  + par.s).^par.sw) +.5;  % q dependent
    v_bias = -par.vb * bn;%par.v * -(sign(data_in.nu_M-.5  - par.s).*abs(data_in.nu_M-.5  - par.s).^par.sv);      % nu dependent
else

    w_bias =  .5;  %par.w *  (sign( data_in.q_M-.5  + par.s).*abs( data_in.q_M-.5  + par.s).^par.sw) +.5;  % q dependent
    v_bias =   0;%par.v * -(sign(data_in.nu_M-.5  - par.s).*abs(data_in.nu_M-.5  - par.s).^par.sv);      % nu dependent
end


% w_bias = par.wq * (qn-.5) + par.wnu * (nun-.5) + .5;  %par.w *  (sign( data_in.q_M-.5  + par.s).*abs( data_in.q_M-.5  + par.s).^par.sw) +.5;  % q dependent
% v_bias = 0; %-par.v * (nun-.5);%par.v * -(sign(data_in.nu_M-.5  - par.s).*abs(data_in.nu_M-.5  - par.s).^par.sv);      % nu dependent

% w_bias = .5;  %par.w *  (sign( data_in.q_M-.5  + par.s).*abs( data_in.q_M-.5  + par.s).^par.sw) +.5;  % q dependent
% v_bias = -par.vnu * (nun-.5) - par.vq * (qn-.5);%par.v * -(sign(data_in.nu_M-.5  - par.s).*abs(data_in.nu_M-.5  - par.s).^par.sv);      % nu dependent


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
end
w_aux =        w_bias + par.kappa_w*data_in.STSE;
w     = 1./(1+exp(-4*(w_aux-.5)));
v     = par.a*(v_bias + par.v_amp*stim_strength + par.kappa_v*data_in.STSE);




p = wfpt_t0noise_vec(reshape(v,1,1,[]),reshape(w,1,1,[]),par.a,par.t0,par.t0_lsig,par.l,o.grid_t,o.k);
% max(sum(p,[1,2]))
% 
% close
% plot(o.grid_t,mean(p(:,2,:),3),'b')
% hold on
% plot(o.grid_t,mean(p(:,1,:),3),'r')
% rtx

% tic
% p_idx = sub2ind(size(p),data_in.rt_norm_idx, data_in.r+1, 1:length(data_in.r));
% toc

nlLH = -sum(log(p(data_in.p_idx)));

if isinf(nlLH)
    nlLH = 10^10;
end

% [data_in.r; data_in.y; p(p_idx)]


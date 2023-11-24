function [LLH distr data_out] = model_LLH_static(par_, par_info, data_in, test_AP, o)


par = parameter_wrap(par_,par_info,'v2s');
par.nu = .5 + (test_AP - par.q);


N_trial = length(data_in.z);
distr = pre_calc_data_generator_static( par, data_in, o);


% % it is taken into account at the end
% % field_names = {'y','z','r'};
% % for ifn = 1:length(field_names)
% %     data_in.(field_names{ifn}) = data_in.(field_names{ifn})(o.interval);
% % end
% % data_in.past_r = data_in.past_r(:,o.interval);

%% noise over trials
auxx = data_in.y;
auxx = (auxx-min(auxx))/(max(auxx)-min(auxx));
auxx = auxx.^par.a;
mu_v = 2*(data_in.z-.5).* (auxx*par.c+par.b);

lP_x_G_tr = -(distr.grid.x-mu_v).^2/2;
lP_x_G_tr = lP_x_G_tr - logsumexp(lP_x_G_tr,1);

%% response probability based on the posterior
% P_r1_G_x = squeeze( 1./(1+exp(-par.beta*(distr.lP_z_G_x(:,:,2) - distr.lP_z_G_x(:,:,1)) + par.kappa*)) );
% P_r1_G_tr = P_r1_G_x' * exp(lP_x_G_tr);


if o.data_gen

    data_out     = data_in;

    past_r = zeros(1,o.n_back);
    w_ST = exp(-par.tau*[1:o.n_back]);
    w_ST = w_ST/sum(w_ST);

    for it=1:N_trial

        %% STSE
        STSE = w_ST*past_r';

%         P_r1_G_x = squeeze( 1./(1+exp(-par.beta*(distr.lP_z_G_x(:,:,2) - distr.lP_z_G_x(:,:,1)) + par.kappa*STSE)) );
        P_r1_G_x = squeeze( 1./(1+exp(-par.beta*(distr.lP_z_G_x(:,:,2) - distr.lP_z_G_x(:,:,1) + par.kappa*STSE))) );
        P_r1_G_tr(it) = P_r1_G_x' * exp(lP_x_G_tr(:,it));


%         %% Lapse
        P_r1_G_tr(it) = P_r1_G_tr(it)*(1-par.lambda)+.5*par.lambda;
% 
        data_out.r(it) = rand < P_r1_G_tr(it);
        data_out.r(it) = P_r1_G_tr(it);

        past_r = [data_out.r(it)*2-1 past_r(1:end-1)];
        % past_r
% 
    end
% 
% %     [data_in.r_v; data_out.r_v]
% 
% 
%     LLH_v = log(P_r1_G_tr(o.interval)).*data_out.r(o.interval) + log(1-P_r1_G_tr(o.interval)).*(1-data_out.r(o.interval));
%     LLH = sum(LLH_v(~isnan(data_out.r(o.interval))));

    LLH=0;


else
    data_out = [];

    %% STSE
    past_r = data_in.past_r*2-1;
    past_r(isnan(past_r)) = 0;
    w_ST = exp(-par.tau*[1:o.n_back]);
    w_ST = w_ST/sum(w_ST);
    STSE = w_ST * past_r;
    % par.tau
    % xxx=w_ST; xxx(1:10)
    % xxx=past_r(:,o.interval); xxx(:,1)'
    % xxx=STSE(o.interval); xxx(1:10)
    % size(STSE)


    % P_r1_G_x_tr = squeeze( 1./(1+exp(-par.beta*(distr.lP_z_G_x(:,:,2) - distr.lP_z_G_x(:,:,1)) + par.kappa*STSE)) );

    % squeeze(distr.lP_z_G_x)
    P_r1_G_x = squeeze( 1./(1+exp(-par.beta*(distr.lP_z_G_x(:,:,2) - distr.lP_z_G_x(:,:,1) + par.kappa*STSE))) );
    P_r1_G_tr = sum(P_r1_G_x .* exp(lP_x_G_tr),1);



    %% Lapse
    P_r1_G_tr = P_r1_G_tr*(1-par.lambda)+.5*par.lambda;
    % P_r1_G_tr(o.interval)
    % size(P_r1_G_tr)


    LLH_v = log(P_r1_G_tr(o.interval)).*data_in.r(o.interval) + log(1-P_r1_G_tr(o.interval)).*(1-data_in.r(o.interval));
    LLH = sum(LLH_v(~isnan(data_in.r(o.interval))));

end 

% LLH=0;


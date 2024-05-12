function distr = pre_calc_data_generator_dynamic(par, distr)

%% grids
distr.grid.x   = linspace(-6,6,distr.N_grid.x)';
distr.grid.y   = linspace(distr.y_min_max(1) , distr.y_min_max(2), distr.N_grid.y);
distr.grid.z   = permute([0 1], [1 ,3,2]);
distr.grid.nu  = permute(linspace(0,1,distr.N_grid.nu+2), [1,3,4,2,5]);   distr.grid.nu = distr.grid.nu(2:end-1);  % q  := AP
distr.grid.q   = permute(linspace(0,1,distr.N_grid.q+2), [1,3,4,5,2]);    distr.grid.q  = distr.grid.q(2:end-1);  % q  := AP


%% P(z|q): prior over z
distr.lP_z_G_q = permute(log([1-distr.grid.q, distr.grid.q]), [1,3,2,4,5]);

%% P(z|q=nu): prior over z
lP_z_G_qEQnu = permute(log([1-distr.grid.nu, distr.grid.nu]), [1,3,2,4,5]);
 
%% P(y): prior over y

distr.lP_y = -2*log(distr.grid.y); % if y=1/gamma
% distr.lP_y = -1.5*log(distr.grid.y); % if y=1/gamma^2 
distr.lP_y = distr.lP_y - logsumexp(distr.lP_y,2);

%% P(x|y,z): noise distribution
mu_v = distr.grid.y.^par.a;
mu_v = (mu_v-min(mu_v))/(max(mu_v)-min(mu_v));
mu_v = (mu_v*par.c+par.b).*(distr.grid.z-.5)*2;

% mu_v = par.c./(1+exp( -par.a*(distr.grid.y-par.b) ));
% mu_v = mu_v.*(distr.grid.z-.5)*2;
distr.mu_v = mu_v;

lP_x_G_y_z = -(distr.grid.x -mu_v).^2/2;
lP_x_G_y_z = lP_x_G_y_z - logsumexp(lP_x_G_y_z, 1);
distr.lP_x_G_y_z = lP_x_G_y_z;
 
%% unbiased P(x|z)
lP_x_G_z_unbiased = logsumexp(lP_x_G_y_z + distr.lP_y, 2);
P_x_G_z_unbiased  = exp(lP_x_G_z_unbiased);
% distr.P_x_G_z_unbiased = P_x_G_z_unbiased;

%% P_{UB}(x) & P_{B}(x|q)
aux_50 = distr.grid.z*0+.5;
P_x_ub_ub     = sum(P_x_G_z_unbiased .* aux_50,3);
P_x_qEQnu_ub  = sum(P_x_G_z_unbiased .* exp(lP_z_G_qEQnu), 3); %% change: exp was missing
% distr.P_x_ub_ub    = P_x_ub_ub;
distr.P_x_qEQnu_ub = P_x_qEQnu_ub;

%% g_{UB}(x)
distr.lP_x_G_z_ratio = lP_x_G_z_unbiased - logsumexp(lP_x_G_z_unbiased, 3);
P_x_G_z_ratio  = exp(distr.lP_x_G_z_ratio);

%% g_{B}(x|nu)
% h = 0.5^.2 - (P_x_G_z_ratio(:,:,2)-0.5).^2;
h = 1 - (2*P_x_G_z_ratio(:,:,2)-1).^2; % change: it has been rewritten
a = (sum(P_x_G_z_ratio(:,:,2) .* P_x_qEQnu_ub, 1) - .5) ./ sum( h .* P_x_qEQnu_ub , 1);
distr.P_x_G_z_ratio_nu = min(1,max(0,( P_x_G_z_ratio - a.*h.*reshape([-1,1],1,1,2) ) ));

%% logP(x|z,nu)
distr.lP_x_G_z_nu = log( P_x_qEQnu_ub .* distr.P_x_G_z_ratio_nu * 2 );
distr.lP_x_G_z_nu = distr.lP_x_G_z_nu - logsumexp(distr.lP_x_G_z_nu,1); % change: add this line to correct for normalization errors

%% logP(x|z,nu,q)
% distr.lP_x_z_G_nu_q = distr.lP_x_G_z_nu + distr.lP_z_G_q;
distr.P_x_z_G_nu_q = exp(distr.lP_x_G_z_nu + distr.lP_z_G_q);

%% logP(x|nu,q)
distr.P_x_G_nu_q = sum(distr.P_x_z_G_nu_q,3);
% distr.lP_x_G_nu_q = logsumexp( distr.lP_x_z_G_nu_q,3);

%% transition probabilites
%% epsilon_q = 1.0/distr["grid"]["q"].size
distr.T_q = -(reshape(distr.grid.q,[],1) - reshape(distr.grid.q,1,[])).^2/2*(10^par.dyn_conc_q);
distr.T_q = distr.T_q-logsumexp(distr.T_q, 1);
distr.T_q = ( exp(distr.T_q) * (par.mix_q) + (1-par.mix_q) * ones(length(distr.grid.q))/length(distr.grid.q) );

distr.T_nu = -(reshape(distr.grid.nu,[],1) - reshape(distr.grid.nu,1,[])).^2/2*(10^par.dyn_conc_nu);
distr.T_nu = distr.T_nu-logsumexp(distr.T_nu, 1);
distr.T_nu = ( exp(distr.T_nu) * (par.mix_nu) + (1-par.mix_nu) * ones(length(distr.grid.nu))/length(distr.grid.nu) );
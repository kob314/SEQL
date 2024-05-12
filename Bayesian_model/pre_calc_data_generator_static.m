function distr = pre_calc_data_generator_static(par, distr)

%% grids
distr.grid.x = linspace(-6,6,distr.N_grid.x)';
distr.grid.y = linspace(distr.y_min_max(1) , distr.y_min_max(2), o.N_grid.y);

%% P(z|q): prior over z
distr.lP_z_G_q = permute(log([1-par.q, par.q]),[1,3,2]);

%% P(y): prior over y
% if ~isfield(par,"sc_y"); par.sc_y = 2; end
% distr.lP_y = -par.sc_y*log(distr.grid.y);
distr.lP_y = -2*log(distr.grid.y);
distr.lP_y = distr.lP_y - logsumexp(distr.lP_y,2);


%% P(x|y,z): noise distribution
auxx = distr.grid.y;
auxx = (auxx-auxx(1))/(auxx(end)-auxx(1)); % scale z between 0 and 1
auxx = auxx.^par.a;
auxx = (auxx*par.c+par.b);
distr.mu_m = auxx.*permute([-1,1],[1,3,2]);


%% log P(x|y,z)
distr.lP_x_G_y_z = -(distr.grid.x-distr.mu_m).^2/2;
distr.lP_x_G_y_z = distr.lP_x_G_y_z - logsumexp(distr.lP_x_G_y_z,1);


%% unbiased P(x|z) = P(x|z;RV=0.5)
lP_x_G_z_nuUB = logsumexp(distr.lP_x_G_y_z + distr.lP_y,2);
P_x_G_z_nuUB  = exp(lP_x_G_z_nuUB);
distr.lP_x_G_z_nuUB = lP_x_G_z_nuUB;

%% P(x|z,RV=nu)
aux_nu = permute([1-par.nu par.nu],[1,3,2]);

%% P(x|RV=0.5,AP=nu) = P(x|RV=nu,AP=0.5) What would be the ground truth P(x) if AP=nu and RV=0.5?
P_z_G_x_nu_qUB  = sum(P_x_G_z_nuUB.*aux_nu,3); % P(x) when q (nu meant to compensate this, that is why aux_nu is the multiplier) is biased but nu is unbiased

%% P(z|x;AP=0.5,RV=0.5)
lP_z_G_x_UB = lP_x_G_z_nuUB - logsumexp(lP_x_G_z_nuUB,3);
P_z_G_x_UB  = exp(lP_z_G_x_UB);
distr.P_z_G_x_UB = P_z_G_x_UB;

%% P(z|x;AP=q,RV=0.5) alternative method
% lP_z_G_x_unbiased = lP_x_G_z_unbiased + distr.lP_z_G_q - logsumexp(lP_x_G_z_unbiased + distr.lP_z_G_q,3);
% P_z_G_x_unbiased  = exp(lP_z_G_x_unbiased);
% distr.P_z_G_x_unbiased = P_z_G_x_unbiased;

%% determining P(x|z,RV)
% h = (1-(2*P_z_G_x_unbiased(:,:,2)-1).^2);
% h = 1-abs(2*P_x_G_z_ratio(:,:,2)-1);
h = 1/4 - (P_z_G_x_UB(:,:,2)-0.5).^2;
a = (sum(P_z_G_x_UB(:,:,2).*P_z_G_x_nu_qUB) - 0.5) / sum(  h .* P_z_G_x_nu_qUB );
a = min(1,max(-1,a));

P_z_G_x_nu_qUB = P_z_G_x_UB - a*h.*permute([-1 1],[1,3,2]);
P_x_G_z_nu = 2*P_z_G_x_nu_qUB .* P_z_G_x_nu_qUB;

% close
% min(P_z_G_x_qUB_nu,[],'all')
% max(P_z_G_x_qUB_nu,[],'all')
% min(aux_nu,[],'all')
% plot(squeeze(P_x_G_z_nu))
% rtx

% P_x_G_z_nu(:,:,2) = P_x_G_z_nu(:,:,2) + min(P_x_G_z_nu(:,:,1),0);
% P_x_G_z_nu(:,:,1) = P_x_G_z_nu(:,:,1) - min(P_x_G_z_nu(:,:,1),0);
% P_x_G_z_nu(:,:,1) = P_x_G_z_nu(:,:,1) + min(P_x_G_z_nu(:,:,2),0);
% P_x_G_z_nu(:,:,2) = P_x_G_z_nu(:,:,2) - min(P_x_G_z_nu(:,:,2),0);
% P_x_G_z_nu        = P_x_G_z_nu./sum(P_x_G_z_nu);


distr.P_x_G_z_nu = P_x_G_z_nu;
distr.lP_x_G_z_nu = log(P_x_G_z_nu);


%% P(z|x): posterior given observation
distr.lP_z_G_x = distr.lP_x_G_z_nu + distr.lP_z_G_q;
distr.lP_z_G_x = distr.lP_z_G_x - logsumexp(distr.lP_z_G_x,3);


 clear all
close all

N_evidence_max = 1000;
dt = 0.003; % after 1000 step it is 3sec
N_trial = 6000;

% analytic

a       = 1;
% v_amp   = 6;
t0      = 0.3;
t0_lsig = .1;
% sc_y    = 1;

% v_subj = reshape([2 4 6],1,1,[]);
% w_subj = .5;
l      = 0;

v= .2;
w=.6;

grid_t = linspace(0,3,N_evidence_max)';
k=21;
p   = wfpt_t0noise_vec(v,w,a,t0,t0_lsig,l,grid_t,k);
p_w = wfpt_vec(v,w,a,grid_t,k);

plot(grid_t,mean(p(:,1,:),3)*N_trial)
hold on
plot(grid_t,mean(p(:,2,:),3)*N_trial)
sum(p,'all')
% plot(grid_t,p_w*N_trial)

% 
% 
% generate data
v_vec = repmat(v,N_evidence_max,N_trial/length(v));
noise = normrnd(0,1,N_evidence_max,N_trial)*20;


cum_evidence = a*w+cumsum(v_vec*dt + noise*dt) -.5*a;
cum_evidence(end,:) = sign(cum_evidence(end,:))*a/2;

[~,RT_idx] = max(abs(cum_evidence)>=a/2);
RT = RT_idx*dt + exp(randn([1,N_trial])*t0_lsig+log(t0));

decision = sign(cum_evidence(sub2ind([N_evidence_max,N_trial],RT_idx,1:N_trial)))';
decision = decision == 1;

% plot(cum_evidence); ylim([-1 1]*.5)
% rt

hs1=histcounts(RT(decision),grid_t);%,'EdgeColor','none');
hs0=histcounts(RT(~decision),grid_t);%,'EdgeColor','none');

%%
plot((grid_t(1:end-1)+grid_t(2:end))/2,smooth(hs1))
plot((grid_t(1:end-1)+grid_t(2:end))/2,smooth(hs0))
% 
% 
%% 
% [g1, g2] = ddm_fpt(mu, bound, delta_t, t_max, ...)
% [g1, g2] = ddm_fpt(v(1),a/2,dt,3);
% 
% plot(grid_t+.3,g1*(grid_t(2)-grid_t(1))*N_trial,'--k')
% plot(grid_t+.3,g2*(grid_t(2)-grid_t(1))*N_trial,'--k')

% [t, b] = ddm_rand_sym(v(1),a/2,dt, N_trial);
% t=t+t0;
% hs1=histcounts(t(b),grid_t);%,'EdgeColor','none');
% hs0=histcounts(t(~b),grid_t);%,'EdgeColor','none');
% plot((grid_t(1:end-1)+grid_t(2:end))/2,smooth(hs1))
% plot((grid_t(1:end-1)+grid_t(2:end))/2,smooth(hs0))

% xlim([.3 1.5])
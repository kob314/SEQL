
clear all
close all

set(0,'DefaultFigureWindowStyle','docked')

t = linspace(0,3,301)';
% t = -t;

v =   .5;
w =  .6;
a =   3; 




k =21;
tic
p = wfpt_vec(v,w,a,t,k);
toc

plot(abs(t),p)


% pe =wfpt_vec(-v,1-w,a,t,k);

% plot(abs(t),p(:,2),'-g')
% hold on
% plot(abs(t),p(:,1),'-r')
% 
% 
norm_check = sum(p,'all')
% max(p)

rt

t0 = .5;
t0_logsig =  .2;
l = 0.1;

tic
p = wfpt_t0noise_vec(v,w,a,t0,t0_logsig,l,t,k);
toc
norm_check = sum(p,'all')

plot(abs(t),p(:,2),'-g')
hold on
plot(abs(t),p(:,1),'-r')

% e =   .01;
% [pc_A,ks_A_v,kl_A_v]=wfpt_Adam(t,v,w,a,t0,e);
% [pe_A,ks_A_v,kl_A_v]=wfpt_Adam(t,-v,1-w,a,t0,e);

% plot(t,pc_A,'g')
% hold on
% plot(t,pe_A,'r')

% %% sanity check
% % ks_v = nan(size(t));
% % for i = 1:length(t)
% %     [p_v(i),ks_v(i),kl_v(i)] = wfpt_NF(t(i),v,a,z*a,e);
% %     [p2_v(i),ks_v(i),kl_v(i)] = wfpt_NF(t(i),-v,a,(1-z)*a,e);
% % %     [p2_v(i),ks_v(i),kl_v(i)] = wfpt_NF(t(i),-v,a,z*a,e);
% % end
% % 
% % plot(t,p_v,'--g')
% % hold on
% % plot(t,p2_v,'--r')
% 
% k = 50;
% 
% tt = (t-t0)/a^2;
% 
% 
% eps_s = 1./sqrt(8*pi*tt) .* exp(-(k-2)^2./(2*tt)); % small t expansion error upper bound
% eps_l = 1./(pi*tt) .* exp(-tt*(k^2*pi^2)/2);   % large t expansion error upper bound
% 
% 
% % use small t expansion when eps_s < eps_l and the small t expansion is
% % good or the large t expansion is bad
% use_small = ( eps_s <= eps_l & tt < (k-1)^2 ) | tt < 1/(pi^2*k^2);
% 
% ks = -k:1:k;
% ps = sum((z+2*ks).*exp(-((z+2*ks).^2)./(2*tt),2) ./sqrt(2*pi*tt.^3);
% 
% kl = 1:k;
% pl = sum(kl.*exp(-(kl.^2)*(pi^2).*tt/2).*sin(kl*pi*z),2)*pi;
% 
% p = pl;
% p(use_small) = ps(use_small);
% p(tt<0) = 0;





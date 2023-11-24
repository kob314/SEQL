function p = wfpt_t0noise_vec(v,w,a,t0,t0_lsig,l,grid_t,k)


%% lognormal non-decision time distribution
p_t0_noise          = lognpdf(grid_t,log(t0),t0_lsig);
p_t0_noise_mass     = 1-lognpdf(grid_t(end),log(t0),t0_lsig);
p_t0_noise          = p_t0_noise/sum(p_t0_noise)*p_t0_noise_mass;
% p_t0_noise(1:end-1) = p_t0_noise(1:end-1)/sum(p_t0_noise(1:end-1))*p_t0_noise_mass;
% p_t0_noise(end)     = 1-p_t0_noise_mass;


p_w = wfpt_vec(v,w,a,grid_t,k);


cond_num = length(v.*w.*a);
p_ = nan(length(grid_t)*2-1,2,cond_num);
for i_cond = [1:cond_num]
    p_(:,1,i_cond) = conv(p_w(:,1,i_cond),p_t0_noise);
    p_(:,2,i_cond) = conv(p_w(:,2,i_cond),p_t0_noise);
end
p_ = [p_(1:(length(grid_t)-1),:,:); sum(p_(length(grid_t):end,:,:),1)];

% hold on
% plot(grid_t,mean(p_w,3),':')
% size(mean(p_,3))
% plot(grid_t,mean(p_,3))
% aux=mean(sum(p_,1),3)
% aux(2)/sum(aux)
% sum(mean(p_,3),'all')



% unif = ones(size(grid_t));
% unif = unif/sum(unif)/2;


% p = p_*(1-l) + unif*l;

% sum(p_t0_noise)
% p = p_*(1-l) + p_t0_noise*l/2;

p_rand = sum(p_,2)/2;
p = p_*(1-l) + p_rand*l;

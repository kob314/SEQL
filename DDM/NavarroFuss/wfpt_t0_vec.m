function p = wfpt_t0_vec(v,w,a,t0,l,grid_t,k)




p_w = wfpt_vec(v,w,a,grid_t,k);

% grid_t_edges = [movmean(grid_t,2); grid_t(end)];

t0_idx = discretize(t0,grid_t);
dt = grid_t(2)-grid_t(1);

w1 = 1-(t0-grid_t(t0_idx))/dt;
w2 = 1-(grid_t(t0_idx+1)-t0)/dt;

p1 = zeros(size(p_w));
p2 = zeros(size(p_w));

p1(t0_idx:end,:) = p_w(1:end-t0_idx+1,:);
p2(t0_idx+1:end,:) = p_w(1:end-t0_idx,:);

p_ = w1*p1 + w2*p2;

% p_ = p_w;

% grid_t(t0_idx)-t0 

% mean(sum(p_w,[1,2]))
% 'x'

% cond_num = length(v.*w.*a);
% p_ = nan(length(grid_t)*2-1,2,cond_num);
% for i_cond = [1:cond_num]
%     p_(:,1,i_cond) = conv(p_w(:,1,i_cond),p_t0_noise);
%     p_(:,2,i_cond) = conv(p_w(:,2,i_cond),p_t0_noise);
% end
% p_ = [p_(1:(length(grid_t)-1),:,:); sum(p_(length(grid_t):end,:,:),1)];

% hold on
% plot(grid_t,mean(p_w,3),':')
% size(mean(p_,3))
% plot(grid_t,mean(p_,3))
% aux=mean(sum(p_,1),3)
% aux(2)/sum(aux)
% sum(mean(p_,3),'all')



unif = ones(size(grid_t));
unif = unif/sum(unif)/2;

p = p_*(1-l) + unif*l;
% p=[];

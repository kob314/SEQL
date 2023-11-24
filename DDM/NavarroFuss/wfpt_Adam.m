function [p ks kl]=wfpt_Adam(t,v,z,a,t0,err)

% original Navarro & Fuss code computes the RT distribution for the error
% trials, with this switch, the reaction time is for the correct trials
v = -v;
z = (1-z)

tt=(t-t0)./(a.^2); % use normalized time

%% calculate number of terms needed for small t and large t
[ks kl] = ks_kl_fun(tt,err)


% if small t is better...
K = ceil(ks);
K = max(K);
K = ceil(K/2);

k = -K:1:K;
ps = sum((z+2*k).*exp(-((z+2*k).^2)/2./tt),2);
ps = ps./sqrt(2*pi*tt.^3);

% else % if large t is better...
K = ceil(kl);
k = 1:K;
pl = sum(k.*exp(-(k.^2)*(pi^2).*tt/2).*sin(k*pi*z),2);
pl = pl*pi;

p = ps;
p(ks>=kl) = pl(ks>=kl);

% convert to f(t|v,a,w)
p=p.*exp(-v.*a.*z -(v.^2).*(t-t0)/2)./(a.^2);

% plot([exp(-v.*a.*z -(v.^2).*(t-t0)/2)./(a.^2) exp(v.*a.*z -(v.^2).*(t-t0)/2)./(a.^2)])
% rt

p(tt<0)=0;

% p = p/sum(p);


% % compute f(tt|0,1,w)
% p=0; %initialize density
% if ks<kl % if small t is better...
%     K=ceil(ks); % round to smallest integer meeting error
%     for k=-floor((K-1)/2):ceil((K-1)/2) % loop over k
%         p=p+(w+2*k)*exp(-((w+2*k)^2)/2/tt); % increment sum
%     end
%     p=p/sqrt(2*pi*tt^3); % add constant term
% 
% else % if large t is better...
%     K=ceil(kl); % round to smallest integer meeting error
%     for k=1:K
%         p=p+k*exp(-(k^2)*(pi^2)*tt/2)*sin(k*pi*w); % increment sum
%     end
%     p=p*pi; % add constant term
% end
% 
% % convert to f(t|v,a,w)
% p=p*exp(-v*a*w -(v^2)*t/2)/(a^2);



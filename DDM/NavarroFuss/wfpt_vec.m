function p = wfpt_vec(v,w,a,grid_t,k)

arguments
    v      (1,1,:) double
    w      (1,1,:) double
    a      (1,1,:) double
    grid_t         double
    k              double
end

cond_num = length(v.*w.*a);

if min(grid_t) == 0
    t = grid_t(2:end);
else
    t = grid_t;
end


tt = t./a^2;

v = [v,  -v];
w = [w, 1-w];


eps_s = 1./sqrt(8*pi*tt) .* exp(-(k-2)^2./(2*tt));  % small t expansion error upper bound
eps_l = 1./(pi*tt) .* exp(-tt*(k^2*pi^2)/2);        % large t expansion error upper bound

use_small = ( eps_s <= eps_l & tt < (k-1)^2 ) | tt < 1/(pi^2*k^2); % when to use small t or large t approx


ks = reshape([-(k-1)/2:1:(k-1)/2],1,1,1,[]);
kl = reshape([1:k],1,1,1,[]);

ps = sum((w+2*ks).*exp(-((w+2*ks).^2)./(2*tt)), 4) ./ sqrt(2*pi*tt.^3);
pl = sum(kl.*exp(-(kl.^2)*(pi^2).*tt/2).*sin(kl*pi.*w), 4)*pi;

p = pl;
p(use_small,:,:) = ps(use_small,:,:);

if min(grid_t) == 0
    p = [zeros(1,2,cond_num); squeeze( p.*exp(-v.*a.*w -(v.^2).*t/2)./(a.^2) )];
else
    p = squeeze( p.*exp(-v.*a.*w -(v.^2).*t/2)./(a.^2) );
end

p = p*(grid_t(2)-grid_t(1));

p = p./max(1,sum(p,[1,2])); % ensure that the integral of probability is not larger then 1, (it can be smaller tough, when RTs ourside of the modelled region has non 0 probability)

% sum(p,'all')
% p = p./sum(p,[1,2]);
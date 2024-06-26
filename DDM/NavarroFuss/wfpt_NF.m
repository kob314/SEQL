function [p ks kl]=wfpt_NF(t,v,a,z,err)

tt=t/(a^2); % use normalized time
w=z/a; % convert to relative start point


% calculate number of terms needed for large t
if pi*tt*err<1 % if error threshold is set low enough
    kl=sqrt(-2*log(pi*tt*err)./(pi^2*tt)); % bound
    kl=max(kl,1/(pi*sqrt(tt))); % ensure boundary conditions met
else % if error threshold set too high
    kl=1/(pi*sqrt(tt)); % set to boundary condition
end

% calculate number of terms needed for small t
if 2*sqrt(2*pi*tt)*err<1 % if error threshold is set low enough
    ks=2+sqrt(-2*tt.*log(2*sqrt(2*pi*tt)*err)); % bound
    ks=max(ks,sqrt(tt)+1); % ensure boundary conditions are met
else % if error threshold was set too high
    ks=2; % minimal kappa for that case
end


% compute f(tt|0,1,w)
p=0; %initialize density
if ks<kl % if small t is better...
    K=ceil(ks); % round to smallest integer meeting error
    for k=-floor((K-1)/2):ceil((K-1)/2) % loop over k
        p=p+(w+2*k)*exp(-((w+2*k)^2)/2/tt); % increment sum
    end
    p=p/sqrt(2*pi*tt^3); % add constant term

else % if large t is better...
    K=ceil(kl); % round to smallest integer meeting error
    for k=1:K
        p=p+k*exp(-(k^2)*(pi^2)*tt/2)*sin(k*pi*w); % increment sum
    end
    p=p*pi; % add constant term
end

% convert to f(t|v,a,w)
p=p*exp(-v*a*w -(v^2)*t/2)/(a^2);

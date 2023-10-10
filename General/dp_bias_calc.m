function [dprime bias mu] = dp_bias_calc(trial_v,resp_v)


idx_nan = isnan(resp_v);
if ~isempty(idx_nan)
    trial_v(idx_nan) = [];
    resp_v(idx_nan)  = [];
end

% if ~isempty(resp_v)
% resp_v
if length(resp_v) > 5
    
    correct_v = 1-abs(trial_v - resp_v);
%     CR = sum(correct_v)/length(correct_v);
%     HR  = (sum( correct_v(logical(trial_v)))     +.5 ) /(sum(trial_v)   +1);    % Hit Rate
%     FAR = (sum( 1-correct_v(logical(1-trial_v))) +.5 ) /(sum(1-trial_v) +1);  % False-Alarm Rate
% d=1;
d=0.0001;
% d=.0000001;
    CR = (sum(correct_v) + d)/(length(correct_v) + 2*d);  
    HR  = (sum( correct_v(logical(trial_v)))      +d ) /(sum(trial_v)  +2*d);   % Hit Rate
    FAR = (sum( 1-correct_v(logical(1-trial_v)))  +d ) /(sum(1-trial_v) +2*d);  % False-Alarm Rate

    

    mu = -norminv(1-CR);  %% theoretical correct ratio curve for an unbiased obs. 
    dprime = norminv(HR)-norminv(FAR);
    bias   = (norminv(HR)+norminv(FAR))/2;
    
else
    
    mu = nan;
    dprime = nan;
    bias   = nan;
    
end

% nHR = norminv(HR)
% nFAR = norminv(FAR)
% bias
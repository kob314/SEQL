function [par_init par_lb par_ub] = parameter_init(par_info,init,lb,ub)

lengths     = cell2mat(struct2cell(par_info));
cum_lengths = [0; cumsum(lengths)];
fn = fields(par_info);

[par_init par_lb par_ub] = deal( nan(1,cum_lengths(end)) );
for n = 1:length(fn)
    par_init((cum_lengths(n)+1):cum_lengths(n+1)) = init(n);
end




% par_init = par_in;
% par_lb   = par_in;
% par_ub   = par_in;
% 
% fn = fields(par_info);
% for n = 1:length(fn)
%     par_in((cum_lengths(n)+1):cum_lengths(n+1)) = par_in.(fn{n});
% end




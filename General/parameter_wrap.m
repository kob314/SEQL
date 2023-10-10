function par_out = parameter_wrap(par_in,par_info,mode)

lengths     = cell2mat(struct2cell(par_info));
cum_lengths = [0; cumsum(lengths)];

fn = fields(par_info);
switch mode
    case 's2v'
        par_out = nan(1,cum_lengths(end));
        for n = 1:length(fn)
            par_out((cum_lengths(n)+1):cum_lengths(n+1)) = par_in.(fn{n});
        end
    case 'v2s'
        par_out =struct();
        for n = 1:length(fn)
            par_out.(fn{n}) = par_in((cum_lengths(n)+1):cum_lengths(n+1));
        end
end

% cum_length = [0 cumsum(par_info.length)];
% switch mode
%     case 's2v'
%         par_out = nan(1,cum_length(end));
%         for n = 1:length(par_info.name)
%             par_out((cum_length(n)+1):cum_length(n+1)) = par_in.(par_info.name{n});
%         end
%     case 'v2s'
%         par_out =struct();
%         for n = 1:length(par_info.name)
%             par_out.(par_info.name{n}) = par_in((cum_length(n)+1):cum_length(n+1));
%         end
% 
% end

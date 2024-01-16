%%% info %%%
% this script was used to create the psychometric bias plots shown in Fig...
% main options
% o.lt_st: showing the total/long-term/short-term biases
% o.plot_num: what experiments to include

clear all
close all

set(0,'DefaultFigureWindowStyle','docked')

%% load data and initilize color settings
include_folders_and_initialize
%%
% load experimental data
load("DATA_EXP3.mat");
data = DATA_EXP;

o.n_back = size(data(1).subj(1).past_r,1);
o.exp_type_v = [1,2,3,4,6,7,14,17]; %experiments to include in the analyzes


o.exp_ord_fit = [1 6 5 4 100 2 3 101 102 103 104 105 106 7 108 109 8];
for k = 1:2  % 1: all; 2: long-term 3: short-term

    o.lt_st = k;

    if o.lt_st == 1
        res=load('data_numpyro_psychom_fit_all_kappa0.mat');
        res2=load('data_numpyro_psychom_fit_2halves_kappa0.mat');
    else
        res=load('data_numpyro_psychom_fit_all.mat');
        res2=load('data_numpyro_psychom_fit_2halves.mat');
    end
    subj_struct = [res.subj_struct(1:2:end); res.subj_struct(2:2:end)]';

    % experiment level bias (mean and std of the hyper-parameter's posterior)
    fit.exp.bias_M  = nan(17,3);
    fit.exp.bias_SD = nan(17,3);

    for ie=[1,2,3,4,6,7,14,17]
        fit.exp.bias_M(ie,:)  = [res.exp.bias(o.exp_ord_fit(ie),1) res2.exp.biasFH(o.exp_ord_fit(ie),1) res2.exp.biasSH(o.exp_ord_fit(ie),1)];
        fit.exp.bias_SD(ie,:) = [res.exp.bias(o.exp_ord_fit(ie),2) res2.exp.biasFH(o.exp_ord_fit(ie),2) res2.exp.biasSH(o.exp_ord_fit(ie),2)];
    end

    % subject level bias (mean and std error of the mean)
    fit_ltst(o.lt_st).subj.bias_M  = nan(17,3);
    fit_ltst(o.lt_st).subj.bias_SE = nan(17,3);
    for i = 1:3
        switch i
            case 1
                aux_bias = res.subj.bias(:,1);
            case 2
                aux_bias = res2.subj.biasFH(:,1);
            case 3
                aux_bias  = res2.subj.biasSH(:,1);
        end
        aux_bias_M = fit.exp.bias_M(subj_struct(:,1)+1,i);

        aux = aux_bias + aux_bias_M;

        fit_ltst(o.lt_st).subj.bias_all(:,i) = aux;

        fit_ltst(o.lt_st).subj.bias_M(:,i)  = accumarray(subj_struct(:,1)+1, aux ,[],@mean,nan);  % mean
        bias_SD = accumarray(subj_struct(:,1)+1, aux ,[],@std,nan);                 % standard deviation
        bias_L  = accumarray(subj_struct(:,1)+1, aux ,[],@length,nan);              % sample size
        fit_ltst(o.lt_st).subj.bias_SE(:,i)  = bias_SD./sqrt(bias_L);                                 % standard error of the mean

    end

end

%% correct for the differences in detectability

% compute the effect of individual objects

% store 1) the exp ID, the 2) rare and 3) frequent object/shape ID and 4) the measured bias - the average bias within the exp group for each trial
ob_pair_m = [];
for ie = o.exp_type_v
    ob_pair   = arrayfun(@(is) [ie data(ie).subj(is).obj_pair], data(ie).incl_subj,'UniformOutput',false);
    ob_pair   = [cell2mat(ob_pair')  fit_ltst(2).subj.bias_all((subj_struct(:,1)+1)==ie,1) - mean( fit_ltst(2).subj.bias_all((subj_struct(:,1)+1)==ie,1) )];
    ob_pair_m = [ob_pair_m; ob_pair];
end

ob_neg    = accumarray( ob_pair_m(:,2), ob_pair_m(:,4),[],@mean,NaN); % bias when the objects are rare 
ob_pos    = accumarray( ob_pair_m(:,3), ob_pair_m(:,4),[],@mean,NaN); % bias when the objects are frequent
ob_neg_SD = accumarray( ob_pair_m(:,2), ob_pair_m(:,4),[],@std,NaN); 
ob_pos_SD = accumarray( ob_pair_m(:,3), ob_pair_m(:,4),[],@std,NaN);
ob_neg_L  = accumarray( ob_pair_m(:,2), ob_pair_m(:,4),[],@length,NaN);
ob_pos_L  = accumarray( ob_pair_m(:,3), ob_pair_m(:,4),[],@length,NaN);

obj_bias.M  = .5*(ob_pos-ob_neg); % the detectability of each object/shape (the shape's contribution to the bias)
obj_bias.SE = .5*sqrt(ob_pos.^2./ob_pos_L+ob_neg.^2./ob_neg_L)

% if o.kappa0 == 1;
%     save('obj_correction_1_2_3_4_6_7_14_17_kappa0.mat','obj_bias')
% else
%     save('obj_correction_1_2_3_4_6_7_14_17.mat','obj_bias')
% end

for ie = o.exp_type_v 
    bias_obj_v((subj_struct(:,1)+1)==ie,1) = arrayfun(@(is) obj_bias.M(data(ie).subj(is).obj_pair(2))-obj_bias.M(data(ie).subj(is).obj_pair(1)), data(ie).incl_subj)';
    % bias_obj_v(:) = 0
end


%%

% fit_ltst(1) = rmfield( fit, "exp");
% fit_ltst(2) = rmfield( fit, "exp");
% fit_ltst(3) = rmfield( fit, "exp");


for k=1:3

    o.lt_st=k; % 1: total; 2: long-term 3: short-term

    fit_ltst(o.lt_st).subj.bias_M  = nan(17,3);
    fit_ltst(o.lt_st).subj.bias_SE = nan(17,3);

    switch o.lt_st

        case {1,2}
            for i = 1:3

                fit_ltst(o.lt_st).subj.bias_all(:,i) = fit_ltst(o.lt_st).subj.bias_all(:,i) - bias_obj_v;
                aux = fit_ltst(o.lt_st).subj.bias_all(:,i);

                fit_ltst(o.lt_st).subj.bias_M(:,i)  = accumarray(subj_struct(:,1)+1, aux ,[],@mean,nan);
                bias_SD = accumarray(subj_struct(:,1)+1, aux ,[],@std,nan);
                bias_L  = accumarray(subj_struct(:,1)+1, aux ,[],@length,nan);
                fit_ltst(o.lt_st).subj.bias_SE(:,i)  = bias_SD./sqrt(bias_L);

            end

        case 3

            for i = 1:3

                switch i
                    case 1
                        aux_tau    = res.global_.tau(1);
                        aux_kappa   = res.global_.kappa(1);
                        incl_      = 81:300;
                    case 2
                        aux_tau    = res2.global_.tau(1);
                        aux_kappa   = res2.global_.kappa(1);
                        incl_      = 81:190;
                    case 3
                        aux_tau    = res2.global_.tau(1);
                        aux_kappa   = res2.global_.kappa(1);
                        incl_      = 191:300;
                end
                w_ST = exp(-aux_tau*[1:30]);
                w_ST = w_ST/sum(w_ST);

                for ie = o.exp_type_v
                    STSE = nan(length(data(ie).incl_subj),1);
                    count=0;
                    for is = data(ie).incl_subj
                        count=count+1;
                        pr = data(ie).subj(is).past_r(:,incl_)*2-1;
                        pr(isnan(pr))=0;
                        STSE(count,:) = mean(aux_kappa*w_ST*pr);
                    end

                    fit_ltst(o.lt_st).subj.bias_all(subj_struct(:,1)+1==ie,i) = STSE;

                    fit_ltst(o.lt_st).subj.bias_M(ie,i)  = nanmean(STSE);
                    fit_ltst(o.lt_st).subj.bias_SE(ie,i) = nanstd(STSE)/sqrt(sum(~isnan(STSE)));
                end

            end
    end

end

save('data_corrected_biases.mat',"fit_ltst","obj_bias","subj_struct")
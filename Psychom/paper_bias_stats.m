%%% info %%%
% this script was used to compute the statistical metrics of psychometric analyzes

clear all
close all

set(0,'DefaultFigureWindowStyle','docked')

%% load data and initilize color settings
include_folders_and_initialize
load('data_corrected_biases.mat',"fit_ltst","obj_bias","subj_struct")

o.exp_type_v = [1,6,7,4,3,2,14,17];

o.lt_st = 2; % 1: total; 2: long-term 3: short-term

%% compute ttests (1: between experiments; 2: between 1st and 2nd half (paired); 3: from 0)
clear ttest_1 ttest_2 ttest_3 stat_
for k=1:3
    for ie = o.exp_type_v

        if k==1
            [h_,p_,ci_,stats_]    = ttest(fit_ltst(o.lt_st).subj.bias_all(subj_struct(:,1)+1==ie,2),fit_ltst(o.lt_st).subj.bias_all(subj_struct(:,1)+1==ie,3));
            cd = abs( meanEffectSize(fit_ltst(o.lt_st).subj.bias_all(subj_struct(:,1)+1==ie,2),fit_ltst(o.lt_st).subj.bias_all(subj_struct(:,1)+1==ie,3),Paired=true,Effect="cohen").Effect );
            ttest_2(ie,:) = deal([stats_.tstat stats_.df p_ cd stats_.df+1]);
        end

        [h_,p_,ci_,stats_]    = ttest(fit_ltst(o.lt_st).subj.bias_all(subj_struct(:,1)+1==ie,k));
        cd = abs( meanEffectSize(fit_ltst(o.lt_st).subj.bias_all( subj_struct(:,1)+1==ie,k),Effect="cohen").Effect );
        ttest_3(ie,k,:) = deal([stats_.tstat stats_.df p_ cd stats_.df+1]);

        for je = o.exp_type_v

            [h_,p_,ci_,stats_] = ttest2(fit_ltst(o.lt_st).subj.bias_all(subj_struct(:,1)+1==ie,k),fit_ltst(o.lt_st).subj.bias_all(subj_struct(:,1)+1==je,k));
            cd = abs( meanEffectSize(fit_ltst(o.lt_st).subj.bias_all(subj_struct(:,1)+1==ie,k),fit_ltst(o.lt_st).subj.bias_all(subj_struct(:,1)+1==je,k),Paired=false,Effect="cohen").Effect );
            ttest_1(ie,je,k,:) = deal([stats_.tstat stats_.df p_ cd]);

        end
    end
end

ttest_between = ttest_1(o.exp_type_v,o.exp_type_v,:,:);
ttest_FHSH    = ttest_2(o.exp_type_v,:);
ttest_from0   = ttest_3(o.exp_type_v,:,:);

%% mean + sem
stat_(:,:,1) = fit_ltst(o.lt_st).subj.bias_M(o.exp_type_v,:);
stat_(:,:,2) = fit_ltst(o.lt_st).subj.bias_SE(o.exp_type_v,:);

%% export results to tables
o.csvwrite_on = 1;

phase_name = {'TOT','FH','SH'};

if o.csvwrite_on
     for k = 1:3
         csvwrite([pwd,'\CSVs\ttest2_t_between_', phase_name{k},'_ltst',num2str(o.lt_st),'.csv'],round(ttest_between(:,:,k,1),2))
         csvwrite([pwd,'\CSVs\ttest2_F_between_', phase_name{k},'_ltst',num2str(o.lt_st),'.csv'],round(ttest_between(:,:,k,2),2))
         csvwrite([pwd,'\CSVs\ttest2_p_between_', phase_name{k},'_ltst',num2str(o.lt_st),'.csv'],round(ttest_between(:,:,k,3),3))
         csvwrite([pwd,'\CSVs\ttest2_Cd_between_',phase_name{k},'_ltst',num2str(o.lt_st),'.csv'],round(ttest_between(:,:,k,4),2))

         csvwrite([pwd,'\CSVs\ttest_from0_',phase_name{k},'_ltst',num2str(o.lt_st),'.csv'],round(ttest_from0(:,k,:),3))
         csvwrite([pwd,'\CSVs\stats_',phase_name{k},'_ltst',num2str(o.lt_st),'.csv'],round(stat_(:,k,:),2))
     end
     csvwrite([pwd,'\CSVs\ttest_FHSH','_ltst',num2str(o.lt_st),'.csv'],round(ttest_FHSH,3))
 end
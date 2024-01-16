clear all
close all
set(0,'DefaultFigureWindowStyle','docked')


curren_dir  = pwd;
idcs   = strfind(curren_dir,'\');
parent_dir = curren_dir(1:idcs(end)-1);
addpath(parent_dir)

data_exps = data_reader;

o.n_back = 30;
o.exp_type_v = [1 6 7 4 3 2 14 15 16 17 18];
o.interval   = 81:300;

DATA_EXP = data_conv(data_exps,o);

save('DATA_EXP3.mat','DATA_EXP')
return
% %%
% ie=18
% X1 = cell2mat( arrayfun(@(is) data_exps.train2.exp_num(ie).subj_num(is).trial_v, 1:length(data_exps.train2.exp_num(ie).subj_num), 'UniformOutput', false)' )
% cats = DATA_EXP(ie).last_jump_type;
% X1a = mean( X1(cats,:) );
% X1b = mean( X1(~cats,:) );
% 
% X2 = mean( cell2mat( arrayfun(@(is) data_exps.test__.exp_num(ie).subj_num(is).trial_v, 1:length(data_exps.test__.exp_num(ie).subj_num), 'UniformOutput', false)' ) );
% 
% close
% plot([1:100],[X1a])
% hold on
% plot([1:100],[X1b])
% plot([101:400],[X2],'Color',[1 1 1]*.5)
% 
% plot([100 100],[0 1],'k')
% plot([1 400]'.*[1 1 1],[1 1]'.*[.5 .65 .75],'k')
% axis([1 400 0 1])
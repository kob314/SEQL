clear all
close all
set(0,'DefaultFigureWindowStyle','docked')

include_folders_and_initialize

data_exps = data_reader;

o.n_back = 30;
o.exp_type_v = [1 6 7 4 3 2 14 15 16 17 18];

DATA_EXP = data_conv(data_exps,o);

save('DATA_EXP.mat','DATA_EXP')

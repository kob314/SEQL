
% fit_handler
clear all
close all
set(0,'DefaultFigureWindowStyle','docked')


%% load in data
include_General
data_ = data_reader;

res  = load("data_numpyro_static_new_subjAP_ysc2.mat"); % load numpyro data with ysc=2 (P(y) ~ y^-ysc; when ysc=2 then it is equal to the Test ground truth)
res2 = load("numpyro_psychom_fit_all.mat");
% res = load("numpyro_static_new (4).mat");
addpath([curren_dir,'\NavarroFuss'])

o.exp_type_v = [1 6 7 4 3 2 14 17];
% o.exp_num_v  = [1 6 5 4 0 2 3];
o.n_back = 30;
data = data_conv(data_,o); %% load in data

o.include  = 81:300;
o.N_grid_rt = 30;                                    % # RT bins
o.bins_rt = linspace(0,3,o.N_grid_rt+1)';            %   RT bins
o.grid_t = (o.bins_rt(1:end-1)+o.bins_rt(2:end))/2;  %   RT grid
% o.k=21;
o.k=15;  % Navarro Fuss Taylor expansion order

%% pooling the data

data_in.e  = cell2mat(arrayfun(@(ie) cell2mat( arrayfun(@(is) data(ie).subj(is).z(o.include)*0+ie, data(ie).incl_subj, 'UniformOutput', false) ), o.exp_type_v, 'UniformOutput', false)); % experimental condition
data_in.s  = cell2mat(arrayfun(@(ie) cell2mat( arrayfun(@(is) data(ie).subj(is).z(o.include)*0+is, data(ie).incl_subj, 'UniformOutput', false) ), o.exp_type_v, 'UniformOutput', false)); % subject id
data_in.z  = cell2mat(arrayfun(@(ie) cell2mat( arrayfun(@(is) data(ie).subj(is).z(o.include), data(ie).incl_subj, 'UniformOutput', false) ), o.exp_type_v, 'UniformOutput', false));      
data_in.r  = cell2mat(arrayfun(@(ie) cell2mat( arrayfun(@(is) data(ie).subj(is).r(o.include), data(ie).incl_subj, 'UniformOutput', false) ), o.exp_type_v, 'UniformOutput', false));
data_in.past_r  = cell2mat(arrayfun(@(ie) cell2mat( arrayfun(@(is) data(ie).subj(is).past_r(:,o.include), data(ie).incl_subj, 'UniformOutput', false) ), o.exp_type_v, 'UniformOutput', false));
y_min = 0.03;
data_in.y       = cell2mat(arrayfun(@(ie) cell2mat( arrayfun(@(is) data(ie).subj(is).y(o.include), data(ie).incl_subj, 'UniformOutput', false) ), o.exp_type_v, 'UniformOutput', false)); % raw y
data_in.y_norm  = cell2mat(arrayfun(@(ie) cell2mat( arrayfun(@(is) (data(ie).subj(is).y(o.include)-min(data(ie).subj(is).y(o.include)))/(max(data(ie).subj(is).y(o.include))-min(data(ie).subj(is).y(o.include)))*(1-y_min) + y_min, data(ie).incl_subj, 'UniformOutput', false) ), o.exp_type_v, 'UniformOutput', false)); % y scaled as in the numpyro code
data_in.y  = data_in.y_norm;
data_in.yz = (data_in.z-.5).*data_in.y*2;
data_in.rt = cell2mat(arrayfun(@(ie) cell2mat( arrayfun(@(is) data(ie).subj(is).rt(o.include), data(ie).incl_subj, 'UniformOutput', false) ), o.exp_type_v, 'UniformOutput', false));
data_in.rt_norm = cell2mat(arrayfun(@(ie) cell2mat( arrayfun(@(is) data(ie).subj(is).rt(o.include)/nanmean(data(ie).subj(is).rt(o.include)), data(ie).incl_subj, 'UniformOutput', false) ), o.exp_type_v, 'UniformOutput', false));

%% throw away out extremely small and large reaction times 
rt_cuts = prctile(data_in.rt_norm,[1 99]);
data_in.r(data_in.rt_norm<max(rt_cuts(1),0)) = nan;
data_in.r(data_in.rt_norm>rt_cuts(2)) = nan;

rt_expcond = accumarray(data_in.e',data_in.rt',[max(o.exp_type_v) 1],@nanmean)';
data_in.rt_norm = data_in.rt_norm.*nanmean(data_in.rt); % the experimental conditions' avg RT is unchanged after all the transformations 


%% compute STSE variable
tau = res.global_.tau(1);%0.17;
weight_ST  = exp(-tau*[1:o.n_back]);
weight_ST  = weight_ST/sum(weight_ST);
past_r = (data_in.past_r-.5)*2;
past_r(isnan(past_r)) = 0;
data_in.STSE = weight_ST*past_r;

%% throw away data where response is nan
idx_excl = isnan(data_in.r);
fn = fields(data_in);
for ifn = 1:length(fn)
    data_in.(fn{ifn}) = data_in.(fn{ifn})(:,~idx_excl);
end


%% discretize time and choice+time distribution
data_in.rt_norm_idx = discretize(data_in.rt_norm,o.bins_rt);
data_in.p_idx = sub2ind([o.N_grid_rt,2,length(data_in.r)],data_in.rt_norm_idx, data_in.r+1, 1:length(data_in.r));



%% subjective AP & RV
load("DATA_EXP4.mat")
subj_AP = cell2mat(arrayfun(@(ie) arrayfun(@(is) DATA_EXP(ie).subj(is).subj_AP, DATA_EXP(ie).incl_subj), o.exp_type_v,'UniformOutput',false))';


count=0;
for ie = o.exp_type_v
    for is = data(ie).incl_subj
        q_M(data_in.e == ie & data_in.s == is) = res.subj.q(ie == (res.subj_struct(:,1)+1) & is == (res.subj_struct(:,2)+1),1);
        nu_M(data_in.e == ie & data_in.s == is) = subj_AP(ie == (res.subj_struct(:,1)+1) & is == (res.subj_struct(:,2)+1)) - q_M(data_in.e == ie & data_in.s == is) + .5;
    end
    shift_(data_in.e==ie) = ~ismember(ie,[2,7]);
end

data_in.q_M    = q_M;
data_in.nu_M   = nu_M;
data_in.shift  = shift_;

%% subjective bias

o.exp_ord_fit = [1 6 5 4 100 2 3 0 0 0 0 0 0 7 0 0 8];
fit.exp.bias_M = nan(7,1);
for ie=o.exp_type_v
   fit.exp.bias_M(ie)  = [res2.exp.bias(o.exp_ord_fit(ie),1)];
end
bias = fit.exp.bias_M(res.subj_struct(:,1)+1,1)+res2.subj.bias(:,1);


for ie = o.exp_type_v
    for is = data(ie).incl_subj
        bias_M(data_in.e == ie & data_in.s == is) = bias(ie == (res.subj_struct(:,1)+1) & is == (res.subj_struct(:,2)+1),1);
    end
end
data_in.bias_M = bias_M;


%%
o.N_grid_t = 30;
o.bins_t = linspace(0,3,o.N_grid_t+1)';
o.grid_t = (o.bins_t(1:end-1)+o.bins_t(2:end))/2;

BIC_=[];
for model = 1:5
    switch model
        case 1
            par_info = struct( ...
                'w'       ,1, ...
                'v'       ,1, ...
                's'       ,1, ...
                'kappa_w' ,1, ...
                'kappa_v' ,1, ...
                'v_amp'   ,1, ...
                'a'       ,1, ...
                't0'      ,1, ...
                't0_lsig' ,1, ...
                'l'       ,1);


            load(['data_res_basedonBayes_complex.mat']);
            par_=par_m(1,:);
            %     par_=[0.4981    3.7545    0.2613    0.0044    0.0748    0.6426   14.1470    1.0901    0.6159    0.1994    0.0622];
            [nlLH(1) pp] = nlLH_DDM_basedonStaticBayes_new(par_,par_info,data_in,o);
            BIC_(model) = length(par_)*log(size(pp,3))+2*nlLH(model);


        case 2

            par_info = struct( ...
                'wb'       ,1, ...
                'vb'       ,1, ...
                's'       ,1, ...
                'kappa_w' ,1, ...
                'kappa_v' ,1, ...
                'v_amp'   ,1, ...
                'a'       ,1, ...
                't0'      ,1, ...
                't0_lsig' ,1, ...
                'l'       ,1);
            load(['data_res_basedonBayes_simple.mat'])
            par_=par_m(1,:);
            [nlLH(2) pp] = nlLH_DDM_basedonStaticBayes_new(par_,par_info,data_in,o);
            BIC_(model) = length(par_)*log(size(pp,3))+2*nlLH(model);


        case 3

            par_info = struct( ...
                'w'       ,1, ...
                'v'       ,1, ...
                's'       ,1, ...
                'v_amp'   ,1, ...
                'a'       ,1, ...
                't0'      ,1, ...
                't0_lsig' ,1, ...
                'l'       ,1);
            load(['data_res_basedonBayes_complex_noPast.mat'])
            par_=par_m(1,:);
            [nlLH(3) pp] = nlLH_DDM_basedonStaticBayes_new(par_,par_info,data_in,o);
            BIC_(model) = length(par_)*log(size(pp,3))+2*nlLH(model);

        case 4

            par_info = struct( ...
                'wb'       ,1, ...
                'vb'       ,1, ...
                's'       ,1, ...
                'v_amp'   ,1, ...
                'a'       ,1, ...
                't0'      ,1, ...
                't0_lsig' ,1, ...
                'l'       ,1);
            load(['data_res_basedonBayes_simple_noPast.mat'])
            par_=par_m(1,:);
            [nlLH(4) pp] = nlLH_DDM_basedonStaticBayes_new(par_,par_info,data_in,o);
            BIC_(model) = length(par_)*log(size(pp,3))+2*nlLH(model);

        case 5

            par_info = struct( ...
                'kappa_w' ,1, ...
                'kappa_v' ,1, ...
                'v_amp'   ,1, ...
                'a'       ,1, ...
                't0'      ,1, ...
                't0_lsig' ,1, ...
                'l'       ,1);
            load(['data_res_basedonBayes_onlyPast.mat'])
            par_=par_m(1,:);
            [nlLH(5) pp] = nlLH_DDM_basedonStaticBayes_new(par_,par_info,data_in,o);
            BIC_(model) = length(par_)*log(size(pp,3))+2*nlLH(model);

    end
end

%%
lLH = -nlLH;
dLH = lLH-min(lLH)

BICp = -BIC_;
dBIC = BICp-min(BICp)

% X = categorical({'complex w STSE','simple w STSE','complex wo STSE','simple wo STSE','only STSE'});
% X = reordercats(X,{'complex w STSE','simple w STSE','complex wo STSE','simple wo STSE','only STSE'});



bar([1:5],dBIC,'EdgeColor','none')
% bar([1:5],dLH,'EdgeColor','none')
text([1:5],[1 1 1 1 1]*1750,string(round(dBIC)),'HorizontalAlignment','center')
axis([.5 5.5 0 1800])
ylabel('\Delta BIC')
xticklabels({'complex w STSE','simple w STSE','complex wo STSE','simple wo STSE','only STSE'})
box off

exportgraphics(gcf,'figures\dBIC.emf')
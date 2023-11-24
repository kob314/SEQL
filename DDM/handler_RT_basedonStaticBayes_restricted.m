% fit_handler

clear all
close all

set(0,'DefaultFigureWindowStyle','docked')

%%
%% load in data

curren_dir  = pwd;
idcs   = strfind(curren_dir,'\');
parent_dir = curren_dir(1:idcs(end)-1);
addpath(parent_dir)

data_ = data_reader;

o.exp_type_v = [1 6 7 4 3 2];
o.exp_num_v  = [1 6 5 4 0 2 3];
o.n_back = 30;
data = data_conv(data_,o);%% load in data

o.include  = 81:300;
o.N_grid_t = 30;
o.bins_t = linspace(0,3,o.N_grid_t+1)';
o.grid_t = (o.bins_t(1:end-1)+o.bins_t(2:end))/2;
% o.k=21;
o.k=15;

%%

data_in.e  = cell2mat(arrayfun(@(ie) cell2mat( arrayfun(@(is) data(ie).subj(is).z(o.include)*0+ie, data(ie).incl_subj, 'UniformOutput', false) ), o.exp_type_v, 'UniformOutput', false));
data_in.s  = cell2mat(arrayfun(@(ie) cell2mat( arrayfun(@(is) data(ie).subj(is).z(o.include)*0+is, data(ie).incl_subj, 'UniformOutput', false) ), o.exp_type_v, 'UniformOutput', false));
data_in.z  = cell2mat(arrayfun(@(ie) cell2mat( arrayfun(@(is) data(ie).subj(is).z(o.include), data(ie).incl_subj, 'UniformOutput', false) ), o.exp_type_v, 'UniformOutput', false));
data_in.r  = cell2mat(arrayfun(@(ie) cell2mat( arrayfun(@(is) data(ie).subj(is).r(o.include), data(ie).incl_subj, 'UniformOutput', false) ), o.exp_type_v, 'UniformOutput', false));
data_in.past_r  = cell2mat(arrayfun(@(ie) cell2mat( arrayfun(@(is) data(ie).subj(is).past_r(:,o.include), data(ie).incl_subj, 'UniformOutput', false) ), o.exp_type_v, 'UniformOutput', false));
data_in.y  = cell2mat(arrayfun(@(ie) cell2mat( arrayfun(@(is) data(ie).subj(is).y(o.include), data(ie).incl_subj, 'UniformOutput', false) ), o.exp_type_v, 'UniformOutput', false));
data_in.y_norm  = cell2mat(arrayfun(@(ie) cell2mat( arrayfun(@(is) data(ie).subj(is).y(o.include)/min(data(ie).subj(is).y(o.include)), data(ie).incl_subj, 'UniformOutput', false) ), o.exp_type_v, 'UniformOutput', false));
% data_in.y_zsc  = cell2mat(arrayfun(@(ie) cell2mat( arrayfun(@(is) (data(ie).subj(is).y - nanmean()), data(ie).incl_subj, 'UniformOutput', false) ), o.exp_type_v, 'UniformOutput', false));
data_in.rt = cell2mat(arrayfun(@(ie) cell2mat( arrayfun(@(is) data(ie).subj(is).rt(o.include), data(ie).incl_subj, 'UniformOutput', false) ), o.exp_type_v, 'UniformOutput', false));
% data_in.rt_zsc = cell2mat(arrayfun(@(ie) cell2mat( arrayfun(@(is) (data(ie).subj(is).rt(o.include) - nanmean(data(ie).subj(is).rt))/nanstd(data(ie).subj(is).rt), data(ie).incl_subj, 'UniformOutput', false) ), o.exp_type_v, 'UniformOutput', false));
% data_in.rt_norm = cell2mat(arrayfun(@(ie) cell2mat( arrayfun(@(is) data(ie).subj(is).rt(o.include) - prctile(data(ie).subj(is).rt,1), data(ie).incl_subj, 'UniformOutput', false) ), o.exp_type_v, 'UniformOutput', false));
data_in.rt_norm = cell2mat(arrayfun(@(ie) cell2mat( arrayfun(@(is) data(ie).subj(is).rt(o.include)/nanmean(data(ie).subj(is).rt(o.include)), data(ie).incl_subj, 'UniformOutput', false) ), o.exp_type_v, 'UniformOutput', false));


% data_in.rt_norm = data_in.rt_norm + nanmean(data_in.rt_norm);
rt_cuts = prctile(data_in.rt_norm,[1 99]);
data_in.r(data_in.rt_norm<max(rt_cuts(1),0)) = nan;
data_in.r(data_in.rt_norm>rt_cuts(2)) = nan;

rt_expcond = accumarray(data_in.e',data_in.rt',[7 1],@nanmean)';
data_in.rt_norm = data_in.rt_norm.*nanmean(data_in.rt);%rt_expcond(data_in.e);
nanmean(data_in.rt)

data_in.y = data_in.y_norm*min(data_in.y);


tau = 0.17;
weight_ST  = exp(-tau*[1:o.n_back]);
weight_ST  = weight_ST/sum(weight_ST);
past_r = (data_in.past_r-.5)*2;
past_r(isnan(past_r)) = 0;
data_in.STSE = weight_ST*past_r;

idx_excl = isnan(data_in.r);
fn = fields(data_in);

for ifn = 1:length(fn)
    data_in.(fn{ifn}) = data_in.(fn{ifn})(:,~idx_excl);
end





% idx_excl = data_in.STSE > 0;
% fn = fields(data_in);
% 
% for ifn = 1:length(fn)
%     data_in.(fn{ifn}) = data_in.(fn{ifn})(:,~idx_excl);
% end
% 
% length(data_in.STSE)


data_in.rt_norm_idx = discretize(data_in.rt_norm,o.bins_t);
data_in.p_idx = sub2ind([o.N_grid_t,2,length(data_in.r)],data_in.rt_norm_idx, data_in.r+1, 1:length(data_in.r));

%
%%
par_info = struct( ...
    'wb'      ,1, ... %'wq'      ,1, ... 'wnu'     ,1, ...
    'vb'      ,1, ...
    'ss'      ,1, ...
    's'       ,1, ...
    'v_amp'   ,1, ...
    'a'       ,1, ...
    't0'      ,1, ...
    't0_lsig' ,1, ...
    'l'       ,1);


% par_info = struct( ...
%     'wb'      ,1, ... %'wq'      ,1, ... 'wnu'     ,1, ...
%     'vb'      ,1, ...
%     'ss'      ,1, ...
%     's'       ,1, ...
%     'kappa_w' ,1, ...
%     'kappa_v' ,1, ...
%     'v_amp'   ,1, ...
%     'a'       ,1, ...
%     't0'      ,1, ...
%     't0_lsig' ,1, ...
%     'l'       ,1);


% par_info = struct( ...
%     'kappa_w' ,1, ...
%     'kappa_v' ,1, ...
%     'v_amp'   ,1, ...
%     'a'       ,1, ...
%     't0'      ,1, ...
%     't0_lsig' ,1, ...
%     'l'       ,1);


% % %      v  w  s tv tw kv kw va  a t0 t0l  l
% l_v = [-3  -3  -.15 0  0 -2 -2  1 .5 .1 .02  0];
% u_v = [ 3   3  .15  5  5  2  2 20  3  1   1 .5];

% %       w   v   ss       s   kw kv va  a t0 t0l  l
% l_v = [ -2 -10    .1   -.15   0  0  1 .5 .1 .02  0];
% u_v = [ 2   10   1.5    .15   2  2 20  3  1 .3 .5];

% % %      wb   vb   ss   s   kw kv va  a t0 t0l  l
% l_v = [ -2  -10   .1   -.5   -2  -2  1 .5 .1 .02  0];
% u_v = [  2   10    3    .5   2  2 20  3  1 .3 .5];

%       wb   vb   ss       s    va  a t0 t0l  l
l_v = [ -2  -10   .1   -.5   1 .5 .1 .02  0];
u_v = [  2   10    3    .5  20  3  1 .3 .5];

%           kw kv va  a t0 t0l  l
% l_v = [   -4  -4  1 .5 .1 .02  0];
% u_v = [    4  4 20  3  1 .3 .5];

% %       w   v   ss       s    va  a t0 t0l  l
% l_v = [ 0   0    .2   -.15     1 .5 .1 .02  0];
% u_v = [ 1   5   1.5    .15    20  3  1 .3 .5];

% %       w   v       s   kv kw va  a t0 t0l  l
% l_v = [ 0   0    -.15   0  0  1 .5 .1 .02  0];
% u_v = [ 3   3     .15   2  2 20  3  1   1 .5];

% %       w   v   sw   sv    s   va  a t0 t0l  l
% l_v = [ 0   0    .1  .1 -.15    1 .5 .1 .02  0];
% u_v = [ 3   3    4    4  .15   20  3  1   1 .5];


%%
ress = load("pymc_static_fit.mat")

% load
resp = load("pymc_psychom_fit.mat")

o.exp_ord_fit = [1 6 5 4 100 2 3];
fit.exp.bias_M = nan(7,1);
for ie=[1,2,3,4,6,7]
   fit.exp.bias_M(ie)  = [resp.exp.bias(o.exp_ord_fit(ie),1)];
end

bias=fit.exp.bias_M(resp.subj_struct(:,1)+1,1)+resp.subj.bias(:,1)

%%
%orig. order 1  2   3   4 5   6  7 
true_AP = [.65 .5 .65 .65 0 .65 .5];
q_M = accumarray(ress.subj_struct(:,1)+1, ress.subj.q(:,1),[],@mean)';
q_L = accumarray(ress.subj_struct(:,1)+1, ress.subj.q(:,1),[],@length)';
q_S = accumarray(ress.subj_struct(:,1)+1, ress.subj.q(:,1),[],@std)'./sqrt(q_L);
nu_M = true_AP-q_M+.5;

xxx.q_M  = q_M;
xxx.nu_M = nu_M;

q_M = nan(1,length(data_in.e));
for ie = o.exp_type_v
    for is = data(ie).incl_subj
        q_M(data_in.e == ie & data_in.s == is) = ress.subj.q(ie == (ress.subj_struct(:,1)+1) & is == (ress.subj_struct(:,2)+1),1);
        nu_M(data_in.e == ie & data_in.s == is) = true_AP(ie)-q_M(data_in.e == ie & data_in.s == is)+.5;
    end
end

data_in.q_M  = q_M;
data_in.nu_M = nu_M;

bias_M = nan(1,length(data_in.e));
for ie = o.exp_type_v
    for is = data(ie).incl_subj
        bias_M(data_in.e == ie & data_in.s == is) = bias(ie == (resp.subj_struct(:,1)+1) & is == (resp.subj_struct(:,2)+1),1);
    end
end

data_in.bias_M = bias_M;



%%
options.UncertaintyHandling = 0;    % Tell BADS that the objective is noisy
options.MaxIter = 200;%150;
options.MaxFunEvals = 1000;

for rep = 1:2

    init_v = l_v;
    for i = 1:length(l_v)
        init_v(i) = rand*(u_v(i)-l_v(i))+l_v(i);
    end

    [par_,fval,exitflag,output] = bads(@(x) nlLH_DDM_basedonStaticBayes_restricted(x,par_info,data_in,o),init_v,l_v,u_v,[],[],[],options);

    par_m(rep,:) = par_
    fval_v(rep)  = fval
end


%%

% save(['res_basedonBayes_simple.mat'],'par_m','fval_v')
save(['res_basedonBayes_simple_noPast.mat'],'par_m','fval_v')
% save(['res_basedonBayes_onlyPast.mat'],'par_m','fval_v')

%%
% init_v = par_m(3,:);
parx = parameter_wrap(par_,par_info,'v2s');

% %%
% close
% o.color = [0.4 0.55 0.25; .45 .25 .5;.9 .35 .4; .2 .3 .6; 0 0 0; .86 .45 .1; 0.93 0.74 0.35];
% 
% subplot(1,2,1)
% hold on
% plot([0 1; .5 .5]',[.5 .5; 0 1]','color',[1 1 1]*0,'LineWidth',1)
% plot([1 0]'+[0 .15],[0 1]'.*[1 1],'color',[1 1 1]*.65,'LineWidth',3)
% % scatter(xxx.q_M(o.exp_type_v) , xxx.nu_M(o.exp_type_v),200,o.color(o.exp_type_v,:),'filled')
% scatter(xxx.q_M(o.exp_type_v) + parx.s, xxx.nu_M(o.exp_type_v) - parx.s,200,o.color(o.exp_type_v,:),'filled')
% xlabel("AP")
% ylabel("ND")
% xticks([.1:.2:.9])
% yticks([.1:.2:.9])
% axis([.25 .75 .25 .75])
% axis([0 1 0 1])
% axis square
% 
% subplot(1,2,2)
% % plot([-1 1; 0 0]',[0 0; -1 1]','color',[1 1 1]*0,'LineWidth',1)
% hold on
% % scatter(parx.w * (xxx.q_M(o.exp_type_v)-.5 + parx.s).^parx.sw , parx.v * (xxx.nu_M(o.exp_type_v)-.5 - parx.s).^parx.sv,200,o.color(o.exp_type_v,:),'filled')
% 
% 
% qnu = [xxx.q_M(o.exp_type_v); xxx.nu_M(o.exp_type_v)];
% qnu = ([1 -1]*qnu).*[1 -1]'/2; % projection
% 
% qnu=qnu*parx.ss; % stretch
% 
% qn = qnu(1,:)+(xxx.q_M(o.exp_type_v)+xxx.nu_M(o.exp_type_v)-1)/2+.5+parx.s; % push back and shift
% nun = qnu(2,:)+(xxx.q_M(o.exp_type_v)+xxx.nu_M(o.exp_type_v)-1)/2+.5-parx.s;
% 
% % scatter(parx.w * (qn-.5) + .5 , parx.v * (nun-.5),200,o.color(o.exp_type_v,:),'filled')
% % scatter(qn ,nun,200,o.color(o.exp_type_v,:),'filled')
% 
% 
% % qn = qnu(1,:)+(xxx.q_M(o.exp_type_v)+xxx.nu_M(o.exp_type_v)-1)/2+.5+parx.s; % push back and shift
% % nun = qnu(2,:)+(xxx.q_M(o.exp_type_v)+xxx.nu_M(o.exp_type_v)-1)/2+.5-parx.s;
% % 
% scatter(parx.w * (qn-.5) + .5 , parx.v * (nun-.5),200,o.color(o.exp_type_v,:),'filled')
% % scatter(qn ,nun,200,o.color(o.exp_type_v,:),'filled')
% 
% % hold on
% % plot([0 1; .5 .5]',[.5 .5; 0 1]','color',[1 1 1]*0,'LineWidth',1)
% % plot([1 0]'+[0 .15],[0 1]'.*[1 1],'color',[1 1 1]*.65,'LineWidth',3)
% % axis([.25 .75 .25 .75])
% % axis([-1 1 -1 1]*.1)
% axis square
% 
% xlabel('w (AP)')
% ylabel('v (ND)')

%%
% load(['res_basedonBayes_additionalScaleParams4.mat']); par_=par_m(1,:);
% par_ = [par_(1:4) 0 0 par_(5:end)]
init_v = par_;
% par_s = parameter_wrap(par_,par_info,'v2s');

% par_current_best = [0.5009    3.7920    0.2594    0.0044    0.0737    0.6457   14.1074    1.0901    0.6159    0.1995    0.0617];
                    %[0.7000    4.4070    0.2105    0.0002    0.0693    0.6786   14.3360    1.0878    0.6178    0.2009    0.0654]
                    %[0.4981    3.7545    0.2613    0.0044    0.0748    0.6426   14.1470    1.0901    0.6159    0.1994    0.0622ú
                    % init_v = [0.3361    2.7522    0.3598    0.0047    0.0833    0.6204   14.0210    1.0929    0.6163    0.1987    0.0577]
% init_v = [0.3361    0    0.3598    0.0047    0.0833    0.6204   14.0210    1.0929    0.6163    0.1987    0.0577]
% init_v = [0    2.7522    0.3598    0.0047    0.0833    0.6204   14.0210    1.0929    0.6163    0.1987    0.0577]
% init_v = [0.3361    2.7522    0.3598    0.0047       14.0210    1.0929    0.6163    0.1987    0.0577]
% % init_v = [0.1858    0.1392    0.2089    0.1444    0.0674    0.1594    7.0173    1.0314    0.5949    0.2095    0.0106];
% init_v(1) = .1; % w
% init_v(2) = .7;  % v
% init_v(3) = .01;  % sw
% init_v(4) = .01;
% % % % % % % % init_v([1 2]) = 0;
% init_v(5) = .05;
% % % % % % % init_v(1) = .2; % w - q dependent
% % % % % % % init_v(2) = .25; % v - nu dependent
% % % % % init_v(end) = .2; 
% % % % % % init_v([6:9]) = [1 1 0 0]; % tau_w tau_v kappa_w kappa_v 
% init_v([5:6]) = [0 0];
% % % init_v([8:9]) = [8 1.5]

% init_v([1 2 6 7 8 9 11 12]) = [.25 5 0 0 10 1.5 .12 .1]
% % init_v
o.N_grid_t = 30;
o.bins_t = linspace(0,3,o.N_grid_t+1)';
o.grid_t = (o.bins_t(1:end-1)+o.bins_t(2:end))/2;
[nlLH pp] = nlLH_DDM_basedonStaticBayes_restricted(init_v,par_info,data_in,o);
nlLH

o.color = [0.4 0.55 0.25; .45 .25 .5;.9 .35 .4; .2 .3 .6; 0 0 0; .86 .45 .1; 0.93 0.74 0.35];

% close
figure
binnum_RT = 3;

count = 0;
for ie = o.exp_type_v
    count = count+1;
   
    
    subplot(6,3,count)
%     % %     plot(o.grid_t,pp1/sum(pp1+pp0),'color',o.color(o.exp_type_v(ie),:),'LineWidth',3)
    %     plot(o.grid_t,smooth(pp1/sum(pp1)),'color',o.color(o.exp_type_v(ie),:),'LineWidth',3)
    %     hold on
    % %     plot(o.grid_t,pp0/sum(pp1+pp0),'LineWidth',3)
    %     plot(o.grid_t,smooth(pp0/sum(pp0)),'LineWidth',3)

    RT       = data_in.rt_norm(data_in.e==ie);
    trial    = data_in.z(data_in.e==ie);
    decision = data_in.r(data_in.e==ie);
    stimstrength = data_in.y(data_in.e==ie) .* (trial-.5)*2;
    correct  = (decision == trial);
% 
    RT_bins = prctile(RT,linspace(0,100,binnum_RT+1));
    RT_cats = discretize(RT,RT_bins);
    RT_cats(isnan(RT_cats)) = binnum_RT+1;
% 
    fc_G_z_M(:,2) = accumarray(RT_cats(trial==1)',correct(trial==1)',[binnum_RT+1 1],@nanmean);
    fc_G_z_M(:,1) = accumarray(RT_cats(trial==0)',correct(trial==0)',[binnum_RT+1 1],@nanmean);
%     fc_G_z_S(:,2) = accumarray(RT_cats(trial==1)',correct(trial==1)',[binnum_RT+1 1],@nanstd);
%     fc_G_z_S(:,1) = accumarray(RT_cats(trial==0)',correct(trial==0)',[binnum_RT+1 1],@nanstd);

    plot(1:3,fc_G_z_M(1:end-1,2)-fc_G_z_M(1:end-1,1),'color',o.color(ie,:),'LineWidth',3)
%     errorbar(1:3,fc_G_z_M(1:end-1,2)-fc_G_z_M(1:end-1,1),mean(fc_G_z_S(1:end-1,:)/sqrt(20),2),'color',o.color(ie,:),'LineWidth',3)
    hold on
    plot([.8 3.2],[0 0],'k')

    axis([.8 3.2 -.3 .42])

    pp_perc = cumsum(mean(sum(pp,2),3));

    [~, pp_perc_idx] = min(abs(pp_perc-[.33 .66]));
    pp_perc_idx = [1 pp_perc_idx length(pp_perc)];


    pp_cf = mean(pp(:,2,data_in.e==ie & data_in.z==1 ),3);
    pp_cr = mean(pp(:,1,data_in.e==ie & data_in.z==0 ),3);
    pp_icf = mean(pp(:,1,data_in.e==ie & data_in.z==1  ),3);
    pp_icr = mean(pp(:,2,data_in.e==ie & data_in.z==0  ),3);

    crf=[];crr=[];
    for i = 1:3
        crf(i) = sum(pp_cf(pp_perc_idx(i):pp_perc_idx(i+1)))/(sum(pp_cf(pp_perc_idx(i):pp_perc_idx(i+1)))+sum(pp_icf(pp_perc_idx(i):pp_perc_idx(i+1))));
        crr(i) = sum(pp_cr(pp_perc_idx(i):pp_perc_idx(i+1)))/(sum(pp_cr(pp_perc_idx(i):pp_perc_idx(i+1)))+sum(pp_icr(pp_perc_idx(i):pp_perc_idx(i+1))));
    end
    

    plot(1:3,crf - crr,'--r','LineWidth',2)


    subplot(6,3,count+6)

    ppe=histcounts(data_in.rt_norm( data_in.e==ie),o.bins_t);
    plot(o.grid_t,ppe/sum(ppe),'color',o.color(ie,:),'LineWidth',3)
    hold on
    %     plot(o.grid_t,ppe/sum(ppe),'color',o.color(o.exp_type_v(ie),:),'LineWidth',3)
    plot(o.grid_t,sum(pp(:,:,data_in.e==ie),[2,3])./sum(pp(:,:,data_in.e==ie),[1,2,3]),'--r','LineWidth',2)


    subplot(6,3,count+12)

    st_edges = [-1 -.3 -.2 -.1  0  .1 .2 .3 1];
    st_cats = discretize(stimstrength,[-1 -.3 -.2 -.1  0  .1 .2 .3 1]);
    psych = accumarray(st_cats',decision',[],@nanmean);

    plot((st_edges(1:end-1)+st_edges(2:end))/2,psych,'color',o.color(ie,:),'LineWidth',3)

    p1 = squeeze(sum(pp(:,2,data_in.e==ie)));
    psych_sim = accumarray(st_cats',p1',[],@nanmean);
    hold on
    plot((st_edges(1:end-1)+st_edges(2:end))/2,psych_sim,'--r','LineWidth',2)

    




end




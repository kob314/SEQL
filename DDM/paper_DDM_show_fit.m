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

%% throw away extremely small and large reaction times 
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

model = 4
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

    case 2

        par_info = struct( ...
            'w'       ,1, ...
            's'       ,1, ...
            'kappa_w' ,1, ...
            'kappa_v' ,1, ...
            'v_amp'   ,1, ...
            'a'       ,1, ...
            't0'      ,1, ...
            't0_lsig' ,1, ...
            'l'       ,1);
        load(['data_res_basedonBayes_simple_AP.mat'])

    case 3

        par_info = struct( ...
            'v'       ,1, ...
            's'       ,1, ...
            'kappa_w' ,1, ...
            'kappa_v' ,1, ...
            'v_amp'   ,1, ...
            'a'       ,1, ...
            't0'      ,1, ...
            't0_lsig' ,1, ...
            'l'       ,1);
        load(['data_res_basedonBayes_simple_RV.mat'])


    case 4

        par_info = struct( ...
            'kappa_w' ,1, ...
            'kappa_v' ,1, ...
            'v_amp'   ,1, ...
            'a'       ,1, ...
            't0'      ,1, ...
            't0_lsig' ,1, ...
            'l'       ,1);
        load(['data_res_basedonBayes_onlyPast.mat'])

end
par_=par_m(1,:);
nlLH(model) = min(fval_v)


par_s = parameter_wrap(par_,par_info,'v2s');
[nlLH pp] = nlLH_DDM_basedonStaticBayes_new(par_,par_info,data_in,o);
nlLH

    


%
o.color = [0.4 0.55 0.25; .45 .25 .5;.9 .35 .4; .2 .3 .6; 0 0 0; .86 .45 .1; 0.93 0.74 0.35];
o.color([14 17],:) = [.6 .3 0; [1 1 1]*.3];
binnum_RT = 3;



t=tiledlayout(3,3,'Units','pixels','InnerPosition',[100 100 800 600]);
% t=tiledlayout(2,4,'Units','pixels','InnerPosition',[100 100 1000 450],'TileSpacing','loose');

% t_v = [3 4 7 8 11 12];
count=0;

fig_order = [1 6 7 4 3 2 14 17];
% fig_order = [1 6 7 14 17];
for ie = fig_order
% for ie = o.exp_type_v
    count=count+1;
%     nexttile(t_v(count))
    nexttile(count)
    hold on
%     [res_(ie).rt_G_rtc_z_sim res_(ie).corr_G_rtc_z_sim] = deal(nan(max(data(ie).incl_subj),n_bins,2));
%     [res_(ie).rt_G_rtc_z_exp res_(ie).corr_G_rtcz_exp]  = deal(nan(max(data(ie).incl_subj),n_bins*2));
    % [rt_G_rtc_z_exp corr_G_rtcz_exp]  = nan(max(data(ie).incl_subj),n_bins);

    [fc_G_z_M sim_fc_G_z_M simw_fc_G_z_M simv_fc_G_z_M sims_fc_G_z_M]= deal(nan(max(data(ie).incl_subj),binnum_RT+1,2));
    for is = data(ie).incl_subj

        RT = data(ie).subj(is).rt(o.include);
        trial    = logical(data(ie).subj(is).z(o.include));
        decision = data(ie).subj(is).r(o.include);
        correct = data(ie).subj(is).correct(o.include);

        RT_bins = prctile(RT,linspace(0,100,binnum_RT+1));
        RT_cats = discretize(RT,RT_bins);
        RT_cats(isnan(RT_cats)) = binnum_RT+1;

        fc_G_z_M(is,:,2) = accumarray(RT_cats(trial)',correct(trial)',[binnum_RT+1 1],@nanmean);
        fc_G_z_M(is,:,1) = accumarray(RT_cats(~trial)',correct(~trial)',[binnum_RT+1 1],@nanmean);

        inc = data_in.e==ie & data_in.s==is;
        p_sim   =  pp(:,:,inc );
%         p_simw   =  ppw(:,:,inc );
%         p_simv   =  ppv(:,:,inc );
%         p_sims   =  pps(:,:,inc );
        pp_perc = cumsum(mean(sum(p_sim,2),3));

        [~, pp_perc_idx] = min(abs(pp_perc-[.33 .66]));
        pp_perc_idx = [1 pp_perc_idx length(pp_perc)];

        pp_cf  = mean(p_sim(:,2,data_in.z(inc)==1 ),3);
        pp_cr  = mean(p_sim(:,1,data_in.z(inc)==0 ),3);
        pp_icf = mean(p_sim(:,1,data_in.z(inc)==1 ),3);
        pp_icr = mean(p_sim(:,2,data_in.z(inc)==0 ),3);

%         ppw_cf  = mean(p_simw(:,2,data_in.z(inc)==1 ),3);
%         ppw_cr  = mean(p_simw(:,1,data_in.z(inc)==0 ),3);
%         ppw_icf = mean(p_simw(:,1,data_in.z(inc)==1 ),3);
%         ppw_icr = mean(p_simw(:,2,data_in.z(inc)==0 ),3);
% 
%         ppv_cf  = mean(p_simv(:,2,data_in.z(inc)==1 ),3);
%         ppv_cr  = mean(p_simv(:,1,data_in.z(inc)==0 ),3);
%         ppv_icf = mean(p_simv(:,1,data_in.z(inc)==1 ),3);
%         ppv_icr = mean(p_simv(:,2,data_in.z(inc)==0 ),3);
% 
%         pps_cf  = mean(p_sims(:,2,data_in.z(inc)==1 ),3);
%         pps_cr  = mean(p_sims(:,1,data_in.z(inc)==0 ),3);
%         pps_icf = mean(p_sims(:,1,data_in.z(inc)==1 ),3);
%         pps_icr = mean(p_sims(:,2,data_in.z(inc)==0 ),3);

        crf=[];crr=[];
        for i = 1:3
            sim_fc_G_z_M(is,i,2) = sum(pp_cf(pp_perc_idx(i):pp_perc_idx(i+1)))/(sum(pp_cf(pp_perc_idx(i):pp_perc_idx(i+1)))+sum(pp_icf(pp_perc_idx(i):pp_perc_idx(i+1))));
            sim_fc_G_z_M(is,i,1) = sum(pp_cr(pp_perc_idx(i):pp_perc_idx(i+1)))/(sum(pp_cr(pp_perc_idx(i):pp_perc_idx(i+1)))+sum(pp_icr(pp_perc_idx(i):pp_perc_idx(i+1))));

%             simw_fc_G_z_M(is,i,2) = sum(ppw_cf(pp_perc_idx(i):pp_perc_idx(i+1)))/(sum(ppw_cf(pp_perc_idx(i):pp_perc_idx(i+1)))+sum(ppw_icf(pp_perc_idx(i):pp_perc_idx(i+1))));
%             simw_fc_G_z_M(is,i,1) = sum(ppw_cr(pp_perc_idx(i):pp_perc_idx(i+1)))/(sum(ppw_cr(pp_perc_idx(i):pp_perc_idx(i+1)))+sum(ppw_icr(pp_perc_idx(i):pp_perc_idx(i+1))));
% 
%             simv_fc_G_z_M(is,i,2) = sum(ppv_cf(pp_perc_idx(i):pp_perc_idx(i+1)))/(sum(ppv_cf(pp_perc_idx(i):pp_perc_idx(i+1)))+sum(ppv_icf(pp_perc_idx(i):pp_perc_idx(i+1))));
%             simv_fc_G_z_M(is,i,1) = sum(ppv_cr(pp_perc_idx(i):pp_perc_idx(i+1)))/(sum(ppv_cr(pp_perc_idx(i):pp_perc_idx(i+1)))+sum(ppv_icr(pp_perc_idx(i):pp_perc_idx(i+1))));
% 
%             sims_fc_G_z_M(is,i,2) = sum(pps_cf(pp_perc_idx(i):pp_perc_idx(i+1)))/(sum(pps_cf(pp_perc_idx(i):pp_perc_idx(i+1)))+sum(pps_icf(pp_perc_idx(i):pp_perc_idx(i+1))));
%             sims_fc_G_z_M(is,i,1) = sum(pps_cr(pp_perc_idx(i):pp_perc_idx(i+1)))/(sum(pps_cr(pp_perc_idx(i):pp_perc_idx(i+1)))+sum(pps_icr(pp_perc_idx(i):pp_perc_idx(i+1))));
        end


    end
    fc_G_z_M = fc_G_z_M(:,1:end-1,:); % remove the values corresponding to nan
    sim_fc_G_z_M = sim_fc_G_z_M(:,1:end-1,:); % remove the values corresponding to nan
    %     simw_fc_G_z_M = simw_fc_G_z_M(:,1:end-1,:);
    %     simv_fc_G_z_M = simv_fc_G_z_M(:,1:end-1,:);
    %     sims_fc_G_z_M = sims_fc_G_z_M(:,1:end-1,:);



    %     plot([0 binnum_RT+1],[.5 .5],'k')
    %     errorbar2(1:binnum_RT, nanmean(fc_G_z_M(:,:,2)), nanstd(fc_G_z_M(:,:,2))/sqrt(length(data(ie).incl_subj)), o.color(ie,:),.3);
    %     plot(1:binnum_RT, nanmean(fc_G_z_M(:,:,2)),'-o','LineWidth',3,'Color',o.color(ie,:),'MarkerFaceColor',o.color(ie,:))
    %     errorbar2(1:binnum_RT, nanmean(fc_G_z_M(:,:,1)), nanstd(fc_G_z_M(:,:,1))/sqrt(length(data(ie).incl_subj)), o.color(ie,:),.3);
    %     ph(2) =plot(1:binnum_RT, nanmean(fc_G_z_M(:,:,1)),'--o','LineWidth',3,'Color',o.color(ie,:),'MarkerFaceColor',o.color(ie,:))
    %
    %     ph(1) =plot(1:binnum_RT, nanmean(sim_fc_G_z_M(:,:,2)),'-o','LineWidth',3,'Color',[1 0 0],'MarkerFaceColor',o.color(ie,:))
    %     plot(1:binnum_RT, nanmean(sim_fc_G_z_M(:,:,1)),'--o','LineWidth',3,'Color',[1 0 0],'MarkerFaceColor',o.color(ie,:))


    plot([0 binnum_RT+1],[0 0],'k','LineWidth',2)

    %     errorbar2(1:binnum_RT, nanmean(fc_G_z_M(:,:,2)-fc_G_z_M(:,:,1)), nanstd(fc_G_z_M(:,:,2)-fc_G_z_M(:,:,1))/sqrt(length(data(ie).incl_subj)), o.color(ie,:),.3);
    %     plot(1:binnum_RT, nanmean(fc_G_z_M(:,:,2)-fc_G_z_M(:,:,1)),'-o','LineWidth',3,'Color',o.color(ie,:),'MarkerFaceColor',o.color(ie,:))

    %     fit_col = [1 0 0];

    %     ph(2) = errorbar(1:binnum_RT, nanmean(fc_G_z_M(:,:,2)-fc_G_z_M(:,:,1)),nanstd(fc_G_z_M(:,:,2)-fc_G_z_M(:,:,1))/sqrt(length(data(ie).incl_subj)),'-o','LineWidth',2,'Color',o.color(ie,:)*0,'MarkerFaceColor',o.color(ie,:)*0,'MarkerEdgeColor',[1 1 1]*0,'MarkerSize',3)

    errorbar2(1:binnum_RT, nanmean(sim_fc_G_z_M(:,:,2)-sim_fc_G_z_M(:,:,1)), nanstd(sim_fc_G_z_M(:,:,2)-sim_fc_G_z_M(:,:,1))/sqrt(length(data(ie).incl_subj)), o.color(ie,:),.3);
    ph(1) = plot(1:binnum_RT, nanmean(sim_fc_G_z_M(:,:,2)-sim_fc_G_z_M(:,:,1)),'-','LineWidth',3,'Color',o.color(ie,:),'MarkerFaceColor',o.color(ie,:));


    %     errorbar(1:binnum_RT, nanmean(fc_G_z_M(:,:,2)-fc_G_z_M(:,:,1)),nanstd(fc_G_z_M(:,:,2)-fc_G_z_M(:,:,1))/sqrt(length(data(ie).incl_subj)),':o','LineWidth',2,'Color',o.color(ie,:)*0,'MarkerFaceColor',o.color(ie,:)*0,'MarkerEdgeColor',[1 1 1]*0,'MarkerSize',4)

    errorbar(1:binnum_RT, nanmean(fc_G_z_M(:,:,2)-fc_G_z_M(:,:,1)),nanstd(fc_G_z_M(:,:,2)-fc_G_z_M(:,:,1))/sqrt(length(data(ie).incl_subj)),':','LineWidth',2,'Color',o.color(ie,:)*0,'capsize',0)
    scatter(1:binnum_RT, nanmean(fc_G_z_M(:,:,2)-fc_G_z_M(:,:,1)),50,1+[0 0 0],'filled','MarkerEdgeColor',[1 1 1]*0,'LineWidth',2.5)
    %     ph(2) = errorbar(1:binnum_RT, nanmean(fc_G_z_M(:,:,2)-fc_G_z_M(:,:,1)),nanstd(fc_G_z_M(:,:,2)-fc_G_z_M(:,:,1))/sqrt(length(data(ie).incl_subj)),'-o','LineWidth',2,'Color',o.color(ie,:)*0,'MarkerFaceColor',o.color(ie,:)*0,'MarkerEdgeColor',[1 1 1]*0,'MarkerSize',3)


    %     plot(1:binnum_RT, nanmean(simw_fc_G_z_M(:,:,2)-simw_fc_G_z_M(:,:,1)),'-','LineWidth',3,'Color',[1 0 0],'MarkerFaceColor',[1 0 0])
    %     plot(1:binnum_RT, nanmean(simv_fc_G_z_M(:,:,2)-simv_fc_G_z_M(:,:,1)),'-','LineWidth',3,'Color',[0 0 1],'MarkerFaceColor',[0 0 1])
    %     plot(1:binnum_RT, nanmean(sims_fc_G_z_M(:,:,2)-sims_fc_G_z_M(:,:,1)),'-','LineWidth',3,'Color',[0 1 0],'MarkerFaceColor',[0 1 0])


    set(gca,'LineWidth',2,'FontSize',22,'FontName','Calibri Light')
    xtickangle(0)

%     title('E')
    box off
%     axis([.9 3.1 .42 .91])
%     axis([.5 3.5 .42 1.1])
%     axis([.5 3.5 -.3 .75])
    axis([.8 3.2 -.3 .4])
    xticks([1 2 3]); xticklabels({'fast','middle','slow'})
%     yticks([.5:.2:.9])
%     axis([.8 3.2 .4 .92])
%     yticks([-.3:.3:.3])

%     if ie == 7
%     legend(ph,{'fit','data'})
%     legend boxoff
%     end


end

%


%

% title(t,'simple model','FontSize',22)
% ylabel(t,'FC bias','FontSize',22)
% xlabel(t,'reaction time category','FontSize',22)


% switch model
%     case 1
%         % title(t,'complex model','FontSize',22)
%         exportgraphics(t,'figures\complex.emf')
%     case 2
%         title(t,'simple model','FontSize',22)
%         exportgraphics(t,'figures\simple2.emf')
%     case 3
%         title(t,'complex model without STSE','FontSize',22)
%         exportgraphics(t,'figures\complex_noPast.emf')
%     case 4
%         title(t,'simple model  without STSE','FontSize',22)
%         exportgraphics(t,'figures\simple_noPast.emf')
%     case 5
%         title(t,'only STSE','FontSize',22)
%         exportgraphics(t,'figures\onlyPast.emf')
% end

% switch model
%     case 1
%         exportgraphics(t,'figures\complex_poster.emf')
% end


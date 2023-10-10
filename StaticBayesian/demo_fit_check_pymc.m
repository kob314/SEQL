%%
clear all
close all

set(0,'DefaultFigureWindowStyle','docked')

%% load in data

curren_dir  = pwd;
idcs   = strfind(curren_dir,'\');
parent_dir = curren_dir(1:idcs(end)-1);
addpath(parent_dir)

data_ = data_reader;

%% load in fit
% load("pymc_static_fit.mat")
% res = load("pymc_static_fit_fixTau.mat")
% res2 = load("numpyro_static (1).mat")
% res2 = load("numpyro_static_old.mat")

% res = load("numpyro_static_new.mat");
res = load("numpyro_static_new (4).mat");
% res = load("numpyro_static_new2.mat");

par.a = res.global_.a(1);
par.b = res.global_.b(1);
par.c = res.global_.c(1);
par.tau   = res.global_.tau(1);
% par.tau   = 0.17;
par.kappa = res.global_.kappa(1);
par.beta  = 100;

% par.a = 0.781;
% par.b = 0.396;
% par.c = 4.0;
% par.kappa = .8
% par.tau   = 0.17;


% % par.a = 1;
% % par.b = 2;
% % par.c = 1;
% par.tau = 1;
% par.kappa = 0;
% par.beta = 100;

par_info = struct( ...
    'a'       ,1, ...
    'b'       ,1, ...
    'c'       ,1, ...
    'tau'     ,1, ...
    'kappa'   ,1, ...
    'beta'    ,1, ...
    'q'       ,1, ...
    'lambda'  ,1);


%%
o.color = [0.4 0.55 0.25; .45 .25 .5;.9 .35 .4; .2 .3 .6; 0 0 0; .86 .45 .1; 0.93 0.74 0.35];
o.color([14 17],:) = [.6 .3 0; [1 1 1]*.3];
o.n_back = 30;
o.exp_type_v = [1 6 7 4 3 2 14 17];% 6 7 4 3 2]
o.interval = 81:300;
data = data_conv(data_,o);

%%
close
o.n_bins_y=4;

data__ = data;
for ie = o.exp_type_v
    for is = data(ie).incl_subj
        data(ie).subj(is).y  = (data(ie).subj(is).y - min(data(ie).subj(is).y(o.interval)))/(max(data(ie).subj(is).y(o.interval)) - min(data(ie).subj(is).y(o.interval)))*(1-0.025) + 0.025;
        data(ie).subj(is).yz = sign(data(ie).subj(is).yz).*data(ie).subj(is).y;
        data(ie).subj(is).y_min_max = [0.025 1];
    end
end
beh_stat = stat_calc(o,data__);

%% grids
% o.N_grid.x = 161;
o.N_grid.x = 100;
o.N_grid.y = 80;
% o.N_grid_nu = 41;
% o.N_grid_q = 41;



LLH_v = 0;
for ie = o.exp_type_v

    data_sim(ie) = data(ie);
    for is = data(ie).incl_subj

        if ~isempty(find((res.subj_struct(:,1)+1==ie) & (res.subj_struct(:,2)+1==is)))
            % [ie, is]
            par.q = res.subj.q((res.subj_struct(:,1)+1==ie) & (res.subj_struct(:,2)+1==is),1);
            par.lambda = res.subj.lambda((res.subj_struct(:,1)+1==ie) & (res.subj_struct(:,2)+1==is),1);

            par_ = parameter_wrap(par,par_info,'s2v');

            o.data_gen = 0;
            [LLH] = model_LLH_static(par_, par_info,  data(ie).subj(is), data(ie).test.AP, o);
            LLH_v(end+1) = LLH;
            % return

            o.data_gen = 1;
            [LLH distr data_sim(ie).subj(is)] = model_LLH_static(par_, par_info,  data(ie).subj(is), data(ie).test.AP, o);
            % LLH_v = LLH_v+LLH;
        else
            data_sim(ie).subj(is).r(:)=nan;
            
        end
    end
end
sum(LLH_v)
% LLH_v

beh_stat_sim = stat_calc(o,data_sim);
% beh_stat_sim

%%
close all
figure(1)
% t=tiledlayout(2,3,'Units','pixels','InnerPosition',[100 100 800 450]);
t=tiledlayout(3,3,'Units','pixels','InnerPosition',[100 100 800 600],'TileSpacing','loose');

% fig_order = [1 6 7 14 4 3 2 17];
fig_order = [1 6 7 4 3 2 14 17];
% fig_order = [1 6 7 14 17];
for ie = fig_order
    ie;
    % ei = ifo;

    nexttile(t)
    hold on
    X = nanmean(beh_stat(ie).yz_M);
    YM = nanmean(beh_stat(ie).r_M);
    YS = nanstd(beh_stat(ie).r_M)/sqrt(length(data(ie).incl_subj));
%     Y = beh_stat(ie).r_M(20,:)
    YM_sim = nanmean(beh_stat_sim(ie).r_M);
    YS_sim = nanstd(beh_stat_sim(ie).r_M)/sqrt(length(data(ie).incl_subj));
%     Y_sim = beh_stat_sim(ie).r_M(20,:)


    plot([-1 1;0 0]',[.5 .5;0 1]','k','LineWidth',1)

%     errorbar2(sign(X).*abs(X).^1,YM,YS, o.color(ie,:),.3);
%     plot(sign(X).*abs(X).^1,YM,'-','LineWidth',3,'Color',o.color(ie,:),'MarkerFaceColor',o.color(ie,:))
%     errorbar(sign(X).*abs(X).^1,YM,YS,'-o','LineWidth',2,'Color',o.color(ie,:)*0,'MarkerFaceColor',o.color(ie,:)*0,'MarkerEdgeColor',o.color(ie,:)*0,'MarkerSize',2) %,'MarkerEdgeColor','none')

    errorbar2(sign(X).*abs(X).^1,YM_sim,YS_sim, o.color(ie,:),.3);
    plot(sign(X).*abs(X).^1,YM_sim,'-','LineWidth',4,'Color',o.color(ie,:))%,'MarkerFaceColor')

    errorbar(sign(X).*abs(X).^1,YM,YS,':','LineWidth',2,'Color',o.color(ie,:)*0,'capsize',0);%,'MarkerFaceColor',o.color(ie,:)*0,'MarkerEdgeColor',o.color(ie,:)*0,'MarkerSize',4) %,'MarkerEdgeColor','none')
    scatter(sign(X).*abs(X).^1,YM,30,[0 0 0]+1,'filled','MarkerEdgeColor',[1 1 1]*0,'LineWidth',2.5)
%     errorbar(sign(X).*abs(X).^1,YM,YS,':','LineWidth',2,'Color',o.color(ie,:)*0,'capsize',0)



%     axis([-.57 .57 0 1]) % 8 bins
    xticks([-.4:.4:.4]), yticks([0:.5:1])
    axis([-.42 .42 0 1]) % 4 bins
    
    set(gca,'LineWidth',2,'FontSize',22,'FontName','Calibri Light')

end

% exportgraphics(t,'figures\fitQual_new.emf')
% exportgraphics(t,'figures\fitQual_new_poster.emf')


%%
% close all
figure(2)
% t=tiledlayout(2,3,'Units','pixels','InnerPosition',[100 100 800 450]);
t=tiledlayout(3,3,'Units','pixels','InnerPosition',[100 100 800 600],'TileSpacing','loose');

% fig_order = [1 6 7 14 4 3 2 17];
fig_order = [1 6 7 4 3 2 14 17];
% fig_order = [1 6 7 14 17];
for ie = fig_order
    ie;
    % ei = ifo;
    % beh_stat(ie).bias_M = beh_stat(ie).bias_M.*sign(beh_stat(ie).dprime_M);
    % beh_stat_sim(ie).bias_M = beh_stat_sim(ie).bias_M.*sign(beh_stat_sim(ie).dprime_M);
    % beh_stat(ie).bias_M = beh_stat(ie).bias_M./abs(beh_stat(ie).dprime_M);
    % beh_stat_sim(ie).bias_M = beh_stat_sim(ie).bias_M./abs(beh_stat_sim(ie).dprime_M);


    nexttile(t)
    hold on
    X = nanmean(beh_stat(ie).y_M);
    YM = nanmean(beh_stat(ie).bias_M);
    YS = nanstd(beh_stat(ie).bias_M)/sqrt(length(data(ie).incl_subj));
%     Y = beh_stat(ie).r_M(20,:2
    YM_sim = nanmean(beh_stat_sim(ie).bias_M);
    YS_sim = nanstd(beh_stat_sim(ie).bias_M)/sqrt(length(data(ie).incl_subj));
%     Y_sim = beh_stat_sim(ie).r_M(20,:)


    plot([-1 1;0 0]',[0 0;0 1]','k','LineWidth',1)

%     errorbar2(sign(X).*abs(X).^1,YM,YS, o.color(ie,:),.3);
%     plot(sign(X).*abs(X).^1,YM,'-','LineWidth',3,'Color',o.color(ie,:),'MarkerFaceColor',o.color(ie,:))
%     errorbar(sign(X).*abs(X).^1,YM,YS,'-o','LineWidth',2,'Color',o.color(ie,:)*0,'MarkerFaceColor',o.color(ie,:)*0,'MarkerEdgeColor',o.color(ie,:)*0,'MarkerSize',2) %,'MarkerEdgeColor','none')

    errorbar2(sign(X).*abs(X).^1,YM_sim,YS_sim, o.color(ie,:),.3);
    plot(sign(X).*abs(X).^1,YM_sim,'-','LineWidth',4,'Color',o.color(ie,:))%,'MarkerFaceColor')

    errorbar(sign(X).*abs(X).^1,YM,YS,':','LineWidth',2,'Color',o.color(ie,:)*0,'capsize',0);%,'MarkerFaceColor',o.color(ie,:)*0,'MarkerEdgeColor',o.color(ie,:)*0,'MarkerSize',4) %,'MarkerEdgeColor','none')
    scatter(sign(X).*abs(X).^1,YM,30,[0 0 0]+1,'filled','MarkerEdgeColor',[1 1 1]*0,'LineWidth',2.5)
%     errorbar(sign(X).*abs(X).^1,YM,YS,':','LineWidth',2,'Color',o.color(ie,:)*0,'capsize',0)



%     axis([-.57 .57 0 1]) % 8 bins
    % xticks([-.4:.4:.4]), yticks([0:.5:1])
    axis([0 .42 -.7 .7]) % 4 bins
    % axis([0 .42 -2 2]) % 4 bins
    
    set(gca,'LineWidth',2,'FontSize',22,'FontName','Calibri Light')

end
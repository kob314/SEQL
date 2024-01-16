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

res=load('data_numpyro_psychom_fit_all.mat');
subj_struct = [res.subj_struct(1:2:end); res.subj_struct(2:2:end)]';

data_ = data_reader;

o.n_back = 30;             % n trials in past to consider
o.exp_type_v = [1 6];      %experiments to include in the analyzes
data = data_conv(data_,o);

o.n_bins_y = 7;
o.interval = [81:300];
beh_stat = stat_calc(o,data); % behavioural data

o.color = [0.4 0.55 0.25; .45 .25 .5;.9 .35 .4; .2 .3 .6; 0 0 0; .86 .45 .1; 0.93 0.74 0.35];
o.color([14 16 17],:) = [1 1 1].*[.5 .3 0]';

%% model

%%% short term serial effects (STSE)%%
n_back=30;
w_ST = exp(-res.global_.tau(1)*[1:n_back]);
w_ST = w_ST./sum(w_ST);

r1_rate = arrayfun(@(ie) mean( arrayfun(@(is) nanmean(data(ie).subj(is).r), data(ie).incl_subj)),o.exp_type_v)';
STSE = res.global_.kappa(1)*sum(w_ST.*r1_rate-w_ST.*(1-r1_rate),2);


%% plot & computation
close
t=tiledlayout(1,1,'Units','pixels','InnerPosition',[130 130 400 500]);
nexttile(t)
hold on

stim_strength =  nan(max(o.exp_type_v),o.n_bins_y*2);
psychom_M     =  nan(max(o.exp_type_v),o.n_bins_y*2);


X = [-.5:.1:.5]; %x axis grid
LW=5;
for count = 1:length(o.exp_type_v)

    if count == 1
        ls = '-';
    else
        ls = '-';
    end

    ie = o.exp_type_v(count); % experiment ID
    stim_strength(ie,:) = nanmean(beh_stat(ie).yz_M(data(ie).incl_subj,:),1);
    psychom_M(ie,:)     = nanmean(beh_stat(ie).r_M(data(ie).incl_subj,:),1);                                    % psychometric data
    psychom_SEM(ie,:)   = nanstd(beh_stat(ie).r_M(data(ie).incl_subj,:),1)/sqrt(length(data(ie).incl_subj));
%     fill([stim_strength(ie,:) fliplr(stim_strength(ie,:))], [psychom_M(ie,:)+psychom_SEM(ie,:) fliplr(psychom_M(ie,:)-psychom_SEM(ie,:))], o.color(ie,:),'EdgeColor','none','FaceAlpha',.2)
%     plot(stim_strength(ie,:), psychom_M(ie,:), 'LineWidth', LW, 'Color', o.color(ie,:))

    subj_idx = (subj_struct(:,1)+1)==ie; % subjects corresponding to the experiment

    % psychometric curve
    Y  = mean( 1./(1+exp(-(res.exp.bias(count,1) + res.subj.bias(subj_idx,1) + res.subj.stim(subj_idx,1).*X + STSE(count) )) ).*(1-res.subj.lambda(subj_idx,1)) + res.subj.lambda(subj_idx,1)*.5);

    ph(count) = plot(X,Y, 'LineWidth', LW, 'Color', o.color(ie,:),'LineStyle',ls);
    sh(count) = scatter(stim_strength(ie,:), psychom_M(ie,:), 80, o.color(ie,:), 'filled','MarkerEdgeColor',o.color(ie,:)*.5);
end

uistack(ph(1),'up',3);
uistack(sh(1),'up',10);

plot([-1 1; 0 0]',[.5 .5; 0 1]','k')

axis([-.5 .5 0 1])
axis([-.6 .5 0 1])
xticks([-1:.5:1])
yticks([0:.5:1])

set(gca,'LineWidth',2,'FontSize',16 )
exportgraphics(t,'figures\psychom_fit2.emf')



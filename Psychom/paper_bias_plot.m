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
load('data_corrected_biases.mat',"fit_ltst","obj_bias","subj_struct")

%% options

% select the type of plot
o.lt_st = 1; % 1: total; 2: long-term 3: short-term
% select what experiments to include
o.plot_num = 5;
switch o.plot_num
    case 1
        o.exp_type_v = [1,6,7];
    case 2
        o.exp_type_v = [1,4,3,2];
    case 3
        o.exp_type_v = [1,14,17];
    case 4
        o.exp_type_v = [1,6,7,4,3,2,14,17];
    case 5
        o.exp_type_v = [1 6]; % Fig 1
end
% select which parts of the test phase are shown 
o.FH_SH = 2; % 0: both, 1: full length, 2: FH + SH
% experiment ordering
o.exp_ord_fit = [1 6 5 4 100 2 3 101 102 103 104 105 106 7 108 109 8]; % correspondance between the original indexing and the indexing used in numpyro

% save figure
o.save_fig_on = 1;



%% plot
close

t = tiledlayout(1,1,'TileSpacing','compact');
t.Units = 'pixels';

t.InnerPosition = [140 120 350 250];
% t.InnerPosition = [140 120 350 150];
if o.plot_num == 1
    t.InnerPosition = [140 120 350 150];
end
nexttile(t)
hold on

FS = 20;    % font size
LW = 14;    % line width

%%%% from a previous version, there is no exp level data analsys anymore %%%%
% o.exp_subj = 2; %1: exp level, 2: subj level
% switch o.exp_subj
%     case 1 % mean: expl level. posterior mean; error bar: posterior std 
%         Y_M = fit.exp.bias_M;
%         Y_S = fit.exp.bias_SD;
%     case 2 % mean: avg. subj level posterior means; error bar: standard error of the mean
%         Y_M = fit.subj.bias_M;
%         Y_S = fit.subj.bias_SE;
% end

Y_M = fit_ltst(o.lt_st).subj.bias_M;
Y_S = fit_ltst(o.lt_st).subj.bias_SE;

N_exp = length(o.exp_type_v);

count = 0;
for ie = o.exp_type_v
    count = count + 1;
    d=-.3;

    switch o.FH_SH
        case 0
        case 1
            d=0;
            w=.4;
            fill(count+[-1 -1 1 1]*w+d,[0 1 1 0]*Y_M(ie,1),o.color(ie,:),'EdgeColor','none')
            errorbar(count+d,Y_M(ie,1),Y_S(ie,2),'.k','LineWidth',3)
%             errorbar(count+d,Y_M(ie,1),Y_S(ie,2),'.','LineWidth',3,'Color',[.8 .8 .8]/es)
        case 2
            d=0;
            w=.4;
            for j = 1:2
                fill(count+[-1 -1 1 1]*w+(N_exp+1)*(j-1),[0 1 1 0]*Y_M(ie,1+j),o.color(ie,:),'EdgeColor',o.color(ie,:)*.7,'linewidth',3)
                errorbar(count+d+(N_exp+1)*(j-1),Y_M(ie,1+j),Y_S(ie,1+j),'.','LineWidth',3,'Color',[0 0 0])
%                 errorbar(count+d+(N_exp+1)*(j-1),Y_M(ie,1+j),Y_S(ie,1+j),'.k','LineWidth',3)
                % errorbar(count+d+(N_exp+1)*(j-1),Y_M(ie,1+j),Y_S(ie,1+j),'.','LineWidth',3,'Color',o.color(ie,:)*.7)
            end
%             plot([0 2*N_exp+2],[0 0],'k','LineWidth',3)
            plot([0 2*N_exp+2],[0 0],'k','LineWidth',3)
            xlim([.3 2*N_exp+1.7])
    end

end

Ymin = min(Y_M(o.exp_type_v,:)-Y_S(o.exp_type_v,:),[],'all');
Ymax = max(Y_M(o.exp_type_v,:)+Y_S(o.exp_type_v,:),[],'all');
ylim([Ymin Ymax]*1.1.*[Ymin<0 Ymax>0]+[0 0]*.9.*[Ymin>0 Ymax<0])

xticklabels([])

ax = gca;
ax.XAxis.Visible = 'off';
set(gca,'FontSize',FS,'FontName','Calibri Light','layer','top','LineWidth',2)

if o.save_fig_on
    switch o.lt_st
        case 1
            switch o.plot_num
                case 1
                    exportgraphics(t,'figures\measured_kappa0_biases_1.emf','ContentType','vector')
                case 2
                    exportgraphics(t,'figures\measured_kappa0_biases_2.emf','ContentType','vector')
                case 3
                    exportgraphics(t,'figures\measured_kappa0_biases_3.emf','ContentType','vector')
                case 5
                    exportgraphics(t,'figures\measured_kappa0_biases_5.emf','ContentType','vector')
            end
        case 2
            switch o.plot_num
                case 1
                    exportgraphics(t,'figures\measured_biases_1.emf','ContentType','vector')
                case 2
                    exportgraphics(t,'figures\measured_biases_2.emf','ContentType','vector')
                case 3
                    exportgraphics(t,'figures\measured_biases_3.emf','ContentType','vector')
                case 5
                    exportgraphics(t,'figures\measured_biases_5.emf','ContentType','vector')
            end
        case 3
            switch o.plot_num
                case 1
                    exportgraphics(t,'figures\measured_STSE_biases_1.emf','ContentType','vector')
                case 2
                    exportgraphics(t,'figures\measured_STSE_biases_2.emf','ContentType','vector')
                case 3
                    exportgraphics(t,'figures\measured_STSE_biases_3.emf','ContentType','vector')
                case 5
                    exportgraphics(t,'figures\measured_STSE_biases_5.emf','ContentType','vector')
            end
    end
end
clear all
close all

set(0,'DefaultFigureWindowStyle','docked')

%% load
load('data_corrected_biases.mat',"fit_ltst","obj_bias","subj_struct")

%% plot
close

t = tiledlayout(1,1,'TileSpacing','compact');
t.Units = 'pixels';
t.InnerPosition = [140 120 350 250];

X  = 1:12;
YM = obj_bias.M;
YS = obj_bias.SE;
[YM,sort_idxs] = sort(YM);

nexttile(t)
hold on

errorbar(1:12,YM,YS(sort_idxs),'.','LineWidth',3)
hold on
plot([0 12],[0 0],'k','LineWidth',3)
xticks([1:11])
yticks([-.2:.1:.2])
xlabel('object ID')
ylabel('object visibility')
axis([0 12 -.2 .2])

set(gca,'FontSize',14,'FontName','Calibri Light','layer','top','LineWidth',1.5)
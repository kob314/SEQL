%% demo_STSE

clear all
close all
set(0,'DefaultFigureWindowStyle','docked')


%% load in data

include_General
data_ = data_reader;

o.exp_type_v = [1 6 7 4 3 2];
o.n_back = 30;
data = data_conv(data_,o);

o.color = [0.4 0.55 0.25; .45 .25 .5;.9 .35 .4; .2 .3 .6; 0 0 0; .86 .45 .1; 0.93 0.74 0.35];

o.n_bins_y = 1;
o.interval = [81:300];
beh_stat = stat_calc(o,data);

close
t=tiledlayout(1,1,'Units','centimeters','InnerPosition',[4 3 18 12]);

LW=3;
nexttile(t)
hold on
STSE = nan(max(o.exp_type_v),o.n_back);
for ie = o.exp_type_v
    STSE(ie,:) = mean(beh_stat(ie).STSE(data(ie).incl_subj,:),1) - .5;
    STSE(ie,:) = movmean(STSE(ie,:),1);
    plot([1:o.n_back], STSE(ie,:), 'LineWidth', LW, 'Color', o.color(ie,:))
end
fill([1 o.n_back o.n_back 1],[-.04 -.04 .08 .08],'w','FaceAlpha',.3)

% fill( [[1:o.n_back] fliplr([1:o.n_back])],[nanmean(STSE)+nanstd(STSE)/sqrt(length(o.exp_type_v)) fliplr(nanmean(STSE)-nanstd(STSE)/sqrt(length(o.exp_type_v)))],'r')
plot([1:o.n_back], nanmean(STSE), 'LineWidth', LW*3, 'Color', [1 .3 1])
errorbar2([1:o.n_back], nanmean(STSE), nanstd(STSE),[1 .3 1],.5)
% /sqrt(length(o.exp_type_v))
plot([1 o.n_back],[0 0],'k','LineWidth',LW)


axis([1 o.n_back -.04 .08])



expEqn = 'c*exp(-tau*x)'
f1 = fit([1:o.n_back]',nanmean(STSE)',expEqn,'Start', [.1 .1])
plot([1:o.n_back]',f1.c*exp(-f1.tau*[1:o.n_back]'),'-g','LineWidth',3)
% plot([1:o.n_back]',.05*exp(-.17*[1:o.n_back]'),'g')

set(gca,'LineWidth',2,'FontSize',26,'FontName','Calibri Light')
box off
xlabel('n back');ylabel('STSE')
exportgraphics(t,'figures\STSE.emf')
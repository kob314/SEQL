clear all
close all

set(0,'DefaultFigureWindowStyle','docked')


o.color = [0.4 0.55 0.25; .45 .25 .5;.9 .35 .4; .2 .3 .6; 0 0 0; .86 .45 .1; 0.93 0.74 0.35];
o.color([14 17],:) = [.6 .35 .2; [1 1 1]*.3];
% o.color([1 6 7],:) = [0.333 0.549 0.341;0.804 0.431 0.235;0.906 0.745 0.373];
% o.color([1 6 7],:) = [0.325 0.478 0.125;0.686 0.333 0.125;0.38 0.596 0.702];
% % o.color([1 6 7],:) = [0.314 0.506 0.016;0.714 0.259 0.004;0.953 0.722 0];
% o.color([1 6 7],:) = [0.4 0.55 0.25; .86 .45 .1; 0.93 0.74 0.35];
% o.color([2 3 4],:) = [.45 .25 .5;.9 .35 .4; .5 .5 .5];
% o.color(4,:) = [.2 .3 .6];
o.based_on_trials_ = 0;  % 1: the statistics depends on the true trial, 0: the statistics depends on the responses

% o.exp_type_v = [2 3 4 1];

o.plot_save_on = 1;
o.total = 1;
o.group_exp = 0; % type of grouping when the dataset is devided into two halfs
o.pred = 1;      % when it is not exp. or sim. data but prediction
o.plot_num = 10;
switch o.plot_num
    case 1
%         o.exp_type_v = [2 3 4 1];
        o.exp_type_v = [1 4 3 2];
    case 2
        o.exp_type_v = [1 14 17];
%     case 2
%         o.exp_type_v = [1 6 7]%;
    case 10
        o.exp_type_v = [1 6 7]%;
end




if o.plot_num == 10
    o.pred_case = 4;
    switch o.pred_case
        case 1
            bias_stat.mean = [];
%             bias_stat.mean([1 6 7],2) = [.45 .25 .28]';   % first half
            bias_stat.mean([1 6 7],2) = [.65 .65 .08]';   % first half
            bias_stat.sem([1 6 7],2) = nan(3,1);
            bias_stat.mean([1 6 7],3) = [.65 .65 .08]';    % second half
            bias_stat.sem([1 6 7],3) = nan(3,1);
        case 2
            
            bias_stat.mean = [];
%             bias_stat.mean([1 6 7],2) = -[.45 .25 .25]';   % first half
            bias_stat.mean([1 6 7],2) = -[.65 .65 .08]';   % first half
            bias_stat.sem([1 6 7],2) = nan(3,1);
            bias_stat.mean([1 6 7],3) = -[.65 .65 .08]';    % second half
            bias_stat.sem([1 6 7],3) = nan(3,1);
        case 3
            bias_stat.mean = [];
%             bias_stat.mean([1 6 7],2) = -[.25 -.25 .25]';   % first half
            bias_stat.mean([1 6 7],2) = -[1 -1 1]*.4';   % first half
            bias_stat.sem([1 6 7],2) = nan(3,1);
            bias_stat.mean([1 6 7],3) = -[1 -1 1]*.4';    % second half
            bias_stat.sem([1 6 7],3) = nan(3,1);
        case 4
            bias_stat.mean = [];
%             bias_stat.mean([1 6 7],2) = [.25 -.25 .25]';   % first half
            bias_stat.mean([1 6 7],2) = [1 -1 1]*.4';   % first half
            bias_stat.sem([1 6 7],2) = nan(3,1);
            bias_stat.mean([1 6 7],3) = [1 -1 1]*.4';    % second half
            bias_stat.sem([1 6 7],3) = nan(3,1);
    end

elseif o.plot_num == 1
    o.pred_case = 2;
    if o.pred_case == 1
        bias_stat.mean = [];
        bias_stat.mean([1 2 3 4],1) = [-.2 .6 .4 .2]';   % first half
        bias_stat.sem( [1 2 3 4],1) = nan(4,1);
    else
        bias_stat.mean = [];
        bias_stat.mean([1 2 3 4],1) = [-.2  .02 -.2 -.2]';   % first half
        bias_stat.sem( [1 2 3 4],1) = nan(4,1);
    end
elseif o.plot_num == 2
    o.pred_case = 1
    if o.pred_case ==1
        bias_stat.mean = [];
        bias_stat.mean([1 14 17],1) = [-.2 -.1 .1]';   % first half
        bias_stat.sem( [1 14 17],1) = nan(3,1);
    else
        bias_stat.mean = [];
        bias_stat.mean([1 14 17],1) = [-.2 -.2 -.2]';   % first half
        bias_stat.sem( [1 14 17],1) = nan(3,1);
    end
end


bias_stat.signif_indiv = [1 1 1]

bias_stat.signif_inc = [];

bias_plot_new3(o, bias_stat)
% axis([.5 4.5 -.3 .7])





group_width = .8;
shift_v = linspace(-1,1,length(o.exp_type_v)*2+1)*group_width/2;
shift_v = shift_v(2:2:end);
bar_half_width = group_width/length(o.exp_type_v)/2;

if o.plot_save_on
    switch o.plot_num
        case 1
            [pwd,'\figures\bias_e1_e4_e3_e2_pred',num2str(o.pred_case),'.emf']
            exportgraphics(gcf,[pwd,'\figures\bias_e1_e4_e3_e2_pred',num2str(o.pred_case),'.emf'])
        case 2
            exportgraphics(gcf,[pwd,'\figures\bias_e1_e14_e17_pred',num2str(o.pred_case),'.emf'])
        case 10
            exportgraphics(gcf,[pwd,'\figures\bias_e1_e6_e7_pred',num2str(o.pred_case),'.emf'])
    end
end

% count = 0;
% for i_e = o.exp_type_v
%     count=count+1;
%     
%     fill([1+shift_v(count)-bar_half_width  1+shift_v(count)+bar_half_width]*[1 1 0 0;0 0 1 1],[0 1 1 0]*bias_stat.mean(i_e,2),o.color(i_e,:),'EdgeColor','none')
%     hold on
%     if ~isnan(bias_stat.sem(i_e,2))
%         errorbar(count-bar_half_width/2,bias_stat.mean(i_e,2),bias_stat.sem(i_e,2),'.','Color','k','LineWidth',4)
%     end
%     
%     fill([2+shift_v(count)-bar_half_width  2+shift_v(count)+bar_half_width]*[1 1 0 0;0 0 1 1],[0 1 1 0]*bias_stat.mean(i_e,3),o.color(i_e,:),'EdgeColor','none')
%     hold on
%     if ~isnan(bias_stat.sem(i_e,2))
%         errorbar(count-bar_half_width/2,bias_stat.mean(i_e,2),bias_stat.sem(i_e,2),'.','Color','k','LineWidth',4)
%     end
%     
% end
% 
% box off
% ax.XAxis.Visible = 'off'
% 
% axis([.5 count+.5 min_Y max_Y])
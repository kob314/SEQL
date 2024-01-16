function bias_plot_new(o, bias_stat)

t = tiledlayout(1,1,'Padding','tight');
t.Units = 'centimeters';
% t.OuterPosition = [1 2 11 10];

if o.plot_num == 1
    t.InnerPosition = [3 2 6 10];
        t.OuterPosition = [1 2 11 10];
elseif o.plot_num == 10
    t.OuterPosition = [1 2 12 6];
else
    t.OuterPosition = [1 2 11 10];
end
nexttile;


exp_name = {['50 ',char(8211),' 65'],['70 ',char(8211),' 50'],['75 ',char(8211),' 65'],['65 ',char(8211),' 65'],'sim',['50  ',char(824),'  65'],['50 ',char(8211),' 65 ', '\' ,' 50'],'50   75   50','50-75 mid shift','50-65Catch','50-65Guess','50-50','50   AS50   50','50 - S50 - 65','50   S75   50','50   V50   65','50   V50   65',['50 ',char(8211),' 65 tilt'],['65 ',char(8211),' 65 tilt']};

count=0;
if o.total
    for i_e = o.exp_type_v
        count=count+1;
        
        bar_half_width = .4;
        [count-bar_half_width count-bar_half_width count+bar_half_width count+bar_half_width]
        fill([count-bar_half_width count-bar_half_width count+bar_half_width count+bar_half_width],[0 1 1 0]*bias_stat.mean(i_e,1),o.color(i_e,:),'EdgeColor',o.color(i_e,:)*.7,'linewidth',3)
        hold on
        %     errorbar(count,bias_stat.mean(i_e),bias_stat.sem(i_e),'.','Color',o.color(i_e,:)*.5,'LineWidth',4)

        if ~isnan(bias_stat.sem(i_e,1))
            errorbar(count,bias_stat.mean(i_e,1),bias_stat.sem(i_e,1),'.','Color','k','LineWidth',4)
        end
        
    end
    
else
    
%     % same experiments next to each other
%     for i_e = o.exp_type_v
%         count=count+1;
%         
%         bar_half_width = .4;
%         fill([count-bar_half_width count-bar_half_width count count],[0 1 1 0]*bias_stat.mean(i_e,2),o.color(i_e,:),'EdgeColor',o.color(i_e,:)*.7,'linewidth',3)
%         hold on
%         if ~isnan(bias_stat.sem(i_e,2))
%             errorbar(count-bar_half_width/2,bias_stat.mean(i_e,2),bias_stat.sem(i_e,2),'.','Color','k','LineWidth',4)
%         end
%         fill([count count count+bar_half_width count+bar_half_width],[0 1 1 0]*bias_stat.mean(i_e,3),.2*o.color(i_e,:)+.8*[1 1 1],'EdgeColor',o.color(i_e,:)*.7,'linewidth',3)
%         
%         if ~isnan(bias_stat.sem(i_e,2))
%             errorbar(count+bar_half_width/2,bias_stat.mean(i_e,3),bias_stat.sem(i_e,3),'.','Color','k','LineWidth',4)
%         end
%         %     errorbar(count,bias_stat.mean(i_e),bias_stat.sem(i_e),'.','Color',o.color(i_e,:)*.5,'LineWidth',4)
%         
%         
%     end
    
    % same periods next to each other
    group_width = .8;
    shift_v = linspace(-1,1,length(o.exp_type_v)*2+1)*group_width/2;
    shift_v = shift_v(2:2:end);
    bar_half_width = group_width/length(o.exp_type_v)/2;
    gap_ = .3;
    for i_e = o.exp_type_v
        count=count+1;
        
        %         fill([1+shift_v(count)-bar_half_width  1+shift_v(count)+bar_half_width]*[1 1 0 0;0 0 1 1],[0 1 1 0]*bias_stat.mean(i_e,2),o.color(i_e,:),'EdgeColor',o.color(i_e,:)*.7,'linewidth',3)
        % bar_half_width
        % shift_v
        %         bias_stat.mean(i_e,2)
%         fill([1+shift_v(count)-bar_half_width  1+shift_v(count)+bar_half_width]*[1 1 0 0;0 0 1 1],[0 1 1 0]*bias_stat.mean(i_e,2),o.color(i_e,:),'EdgeColor','none')
        hold on
        fill([1+shift_v(count)-bar_half_width  1+shift_v(count)+bar_half_width]*[1 1 0 0;0 0 1 1],[0 1 1 0]*bias_stat.mean(i_e,2),o.color(i_e,:),'EdgeColor',o.color(i_e,:)*.7,'linewidth',3)
        if ~isnan(bias_stat.sem(i_e,2))
%             errorbar(count-bar_half_width/2,bias_stat.mean(i_e,2),bias_stat.sem(i_e,2),'.','Color','k','LineWidth',4)
            errorbar(1+shift_v(count),bias_stat.mean(i_e,2),bias_stat.sem(i_e,2),'.','Color','k','LineWidth',4)
        end
        
%         fill([(1+group_width+gap_)+shift_v(count)-bar_half_width  (1+group_width+gap_)+shift_v(count)+bar_half_width]*[1 1 0 0;0 0 1 1],[0 1 1 0]*bias_stat.mean(i_e,3),o.color(i_e,:),'EdgeColor','none')
        fill([(1+group_width+gap_)+shift_v(count)-bar_half_width  (1+group_width+gap_)+shift_v(count)+bar_half_width]*[1 1 0 0;0 0 1 1],[0 1 1 0]*bias_stat.mean(i_e,3),o.color(i_e,:),'EdgeColor',o.color(i_e,:)*.7,'linewidth',3)
        hold on
        if ~isnan(bias_stat.sem(i_e,3))
%             errorbar(count-bar_half_width/2,bias_stat.mean(i_e,3),bias_stat.sem(i_e,3),'.','Color','k','LineWidth',4)
            errorbar((1+group_width+gap_)+shift_v(count),bias_stat.mean(i_e,3),bias_stat.sem(i_e,3),'.','Color','k','LineWidth',4)
        end

    end
    
end



if o.total
    bias_M = bias_stat.mean(:,1);
    bias_S = bias_stat.sem(:,1);
    bias_S(isnan(bias_S)) = 0;
else
    bias_M = bias_stat.mean(:,[2,3]);
    bias_S = bias_stat.sem(:,[2,3]);
    bias_M = bias_M(:);
    bias_S = bias_S(:);
    bias_S(isnan(bias_S)) = 0;
end
if min(bias_M-bias_S) < 0
    min_Y = min(bias_M-bias_S)*1.1;
else
    min_Y = min(bias_M-bias_S)*.9;
end
max_y = max(bias_M+bias_S);
if max(bias_M+bias_S) > 0
    max_Y = max_y*1.1;
else
    max_Y = max_y*.9;
end
d_ = max_Y - max_y;

d_ = .1*(max_Y-min_Y);
% d_ = d_*2;
LWS = 4;
signif_inc = bias_stat.signif_inc;
Y_maxpos = max_Y;
% 
% if o.total
%     for i_signif = 1:size(signif_inc,1)
%         sf = bias_stat.signif_total(o.exp_type_v(signif_inc(i_signif,1)),o.exp_type_v(signif_inc(i_signif,2)));
%         
%         
%         %% version without n.s.
%         if sf < .05
% %             w_ = .2;
%             w_ = .1;
%             Y_pos = max_y + d_*signif_inc(i_signif,3);
%             plot([[1 1]*(signif_inc(i_signif,1)+w_) [1 1]*(signif_inc(i_signif,2)-w_)],[Y_pos-d_/2 [1 1]*Y_pos Y_pos-d_/2],'k','LineWidth',LWS)
%             Y_maxpos = max(Y_maxpos,Y_pos)
%             if sf < .01
%                 text(mean(signif_inc(i_signif,1:2)),Y_pos+d_/4,'\bf**','FontSize',30,'HorizontalAlignment','center')
%             else
%                 text(mean(signif_inc(i_signif,1:2)),Y_pos+d_/4,'\bf*','FontSize',30,'HorizontalAlignment','center')
%             end
%         end
%         
%         %     w_ = .2;
%         %     Y_pos = max_y + d_*signif_inc(i_signif,3);
%         %     plot([[1 1]*(signif_inc(i_signif,1)+w_) [1 1]*(signif_inc(i_signif,2)-w_)],[Y_pos-d_/2 [1 1]*Y_pos Y_pos-d_/2],'k','LineWidth',LWS)
%         %     Y_maxpos = max(Y_maxpos,Y_pos)
%         %     if sf < .05
%         %         if sf < .01
%         %             text(mean(signif_inc(i_signif,1:2)),Y_pos+d_/4,'\bf**','FontSize',30,'HorizontalAlignment','center')
%         %         else
%         %             text(mean(signif_inc(i_signif,1:2)),Y_pos+d_/4,'\bf*','FontSize',30,'HorizontalAlignment','center')
%         %         end
%         %     else
%         %         text(mean(signif_inc(i_signif,1:2)),Y_pos+d_/2,'\bfn.s.','FontSize',20,'HorizontalAlignment','center')
%         %     end
%         
%         %% version with ns
%     end
% else
%     for i_signif = 1:size(signif_inc,1)
%         sf = bias_stat.signif_halfs(o.exp_type_v(signif_inc(i_signif,1)));
%         w_=.05;
%         
%         Y_pos = max_y + d_*signif_inc(i_signif,2);
%         signif_inc(i_signif,1)-bar_half_width/2
%         plot([[1 1]*(signif_inc(i_signif,1)-bar_half_width/2+w_) [1 1]*(signif_inc(i_signif,1)+bar_half_width/2-w_)],[Y_pos-d_/2 [1 1]*Y_pos Y_pos-d_/2],'k','LineWidth',LWS)
%         Y_maxpos = max(Y_maxpos,Y_pos);
%         if sf < .05
%             if sf < .01
%                 text(signif_inc(i_signif,1),Y_pos+d_/4,'\bf**','FontSize',30,'HorizontalAlignment','center')
%             else
%                 text(signif_inc(i_signif,1),Y_pos+d_/4,'\bf*','FontSize',30,'HorizontalAlignment','center')
%             end
%         else
%             text(signif_inc(i_signif,1),Y_pos+d_/2,'\bfn.s.','FontSize',20,'HorizontalAlignment','center')
%         end
%         
%         % version with ns
%     end
% end




max_Y = Y_maxpos;

if o.total | o.group_exp
    axis([.5 count+.6 min_Y max_Y])
    xticks([1:count])
    plot([.5 count+.6],[0 0],'k','LineWidth',3)
else
    axis([.5 (1+group_width+gap_)+.5  min_Y max_Y])
    xticks([1 2])
    plot([.5 (1+group_width+gap_)+.5],[0 0],'k','LineWidth',3)
end
% xticklabels(exp_name(o.exp_type_v))
xticklabels([])
% yticks([-1:.1:1])
% yticks([-1.5:.3:1.5])
if o.pred
    yticks([])
else
%     yticks([-1.2:.4:1.2])
    yticks([-.5:.1:.5])
end

box off



ax = gca;
if o.plot_num == 1

%    axis([.5 4.5 -.5 1.2])
%    axis([.5 4.5 -.5 .9])
   axis([.5 4.5 -.1 .22])
    
   bias_stat.signif_indiv
   count = 0;
% %    for i_e = o.exp_type_v
% %        count = count+1;
% %        
% %        o.color_ND = [0.8320  0.1563  0.0195];
% %        
% %        y_coord = -.55;
% %        if bias_stat.signif_indiv(i_e) < 0.01
% %            if i_e ~= 1
% %                text(count,y_coord,'\bf**','FontSize',30,'HorizontalAlignment','center','VerticalAlignment','middle')
% %            else
% %                text(count,y_coord,'\bf**','FontSize',30,'HorizontalAlignment','center','VerticalAlignment','middle','Color',o.color_ND)
% %            end
% %        elseif bias_stat.signif_indiv(i_e) < 0.05
% %            text(count,y_coord,'\bf*','FontSize',30,'HorizontalAlignment','center','VerticalAlignment','middle')
% %        else
% %            text(count,y_coord+.06,'n.s.','FontSize',25,'HorizontalAlignment','center','VerticalAlignment','middle','Color',o.color_ND)
% %        end 
% %    end
% % 
% %    
% %    %    axis([.5 4.5 -.4 .2])
% %    % plot([[1.1 1.1 1.9 1.9]; [2.1 2.1 2.9 2.9]; [1.1 1.1 2.9 2.9]]',[[0 .07 .07  0]+1.3; [0 .07 .07  0]+0; [0 .07 .07  0]+1.7]','k','LineWidth',LWS)
% %    
% %    
% %    %%exp
% % % %    plot([[1.1 1.1 1.9 1.9]; [2.1 2.1 2.9 2.9]; [2.1 2.1 3.9 3.9]]',[[0 .07 .07  0]+1.3; [0 .07 .07  0]+1.3; [0 .07 .07  0]+1]','k','LineWidth',LWS)
% % %    plot([[2.1 2.1 2.9 2.9]; [2.1 2.1 3.9 3.9]]',[[0 .07 .07  0]+1.3; [0 .07 .07  0]+1]','k','LineWidth',LWS)
% % %    plot([[1.1 1.1 1.9 1.9]]',[[0 .07 .07  0]+1.3]','LineWidth',LWS,'Color',o.color_ND)
% % %    
% % %    text(1.5,1.42,'\bf*','FontSize',25,'HorizontalAlignment','center','VerticalAlignment','middle','Color',o.color_ND)
% % %    text(2.5,1.42,'\bf*','FontSize',25,'HorizontalAlignment','center','VerticalAlignment','middle')
% % %    text(3,1.12,'\bf**','FontSize',25,'HorizontalAlignment','center','VerticalAlignment','middle')
% % %    
% %    
% % 
% %    %%sim
% %    plot([[2.1 2.1 2.9 2.9]; [3.1 3.1 3.9 3.9]]',[[0 .07 .07  0]+1.3; [0 .07 .07  0]+1.3]','k','LineWidth',LWS)
% %    plot([[1.1 1.1 1.9 1.9]]',[[0 .07 .07  0]+1.3]','LineWidth',LWS,'Color',o.color_ND)
% %    
% %    text(1.5,1.42,'\bf**','FontSize',25,'HorizontalAlignment','center','VerticalAlignment','middle','Color',o.color_ND)
% % %    text(2.5,1.42,'\bf**','FontSize',25,'HorizontalAlignment','center','VerticalAlignment','middle')
% % %    text(3.5,1.42,'\bf**','FontSize',25,'HorizontalAlignment','center','VerticalAlignment','middle')
% % 
% % 
   xticklabels(exp_name(o.exp_type_v))
   ax.XAxis.Visible = 'off';
   
elseif o.plot_num == 2
   

    
    
%     axis([.5 3.5 -.4 .2])
%     plot([[1.1 1.1 1.9 1.9]; [2.1 2.1 2.9 2.9]; [1.1 1.1 2.9 2.9]]',[[0 .05 .05  0]+.22; [0 .05 .05  0]+.22; [0 .05 .05  0]+.34]','k','LineWidth',LWS)
% %     text(1.5,.3,'\bf**','FontSize',25,'HorizontalAlignment','center','VerticalAlignment','middle')
% %     text(2.5,.3,'\bf*','FontSize',25,'HorizontalAlignment','center','VerticalAlignment','middle')
% %     text(2,.46,'n.s.','FontSize',25,'HorizontalAlignment','center','VerticalAlignment','middle')
%     
%     text(1.5,.3,'\bf*','FontSize',25,'HorizontalAlignment','center','VerticalAlignment','middle')
%     text(2.5,.3,'\bf*','FontSize',25,'HorizontalAlignment','center','VerticalAlignment','middle')
%     text(2,.46,'n.s.','FontSize',25,'HorizontalAlignment','center','VerticalAlignment','middle')
%     
%     xticklabels(exp_name(o.exp_type_v))
    ax.XAxis.Visible = 'off';

elseif o.plot_num == 4
    ax.XAxis.Visible = 'off';
elseif o.plot_num == 9
    
    %     text(1,.05,'\bf**','FontSize',25,'HorizontalAlignment','center','VerticalAlignment','middle')
    %     text(2,.05,'\bf*','FontSize',25,'HorizontalAlignment','center','VerticalAlignment','middle')
    %     text(count,y_coord+.06,'n.s.','FontSize',25,'HorizontalAlignment','center','VerticalAlignment','middle')
    
    count = 0;
    for i_e = o.exp_type_v
        count = count+1;
        y_coord = -.5;%0.05;
        if bias_stat.signif_indiv(i_e) < 0.01
            text(count,y_coord,'\bf**','FontSize',30,'HorizontalAlignment','center','VerticalAlignment','middle')
        elseif bias_stat.signif_indiv(i_e) < 0.05
            text(count,y_coord,'\bf*','FontSize',30,'HorizontalAlignment','center','VerticalAlignment','middle')
        else
            text(count,y_coord+.06,'n.s.','FontSize',25,'HorizontalAlignment','center','VerticalAlignment','middle')
        end
    end
    
    ax.XAxis.Visible = 'off'
    plot([1 1 2 2]',[[0 .05 .05  0]+.05]','k','LineWidth',LWS)
    text(1.5,.15,'n.s.','FontSize',25,'HorizontalAlignment','center','VerticalAlignment','middle')
    
elseif o.plot_num == 10
    
    %     text(1,.05,'\bf**','FontSize',25,'HorizontalAlignment','center','VerticalAlignment','middle')
    %     text(2,.05,'\bf*','FontSize',25,'HorizontalAlignment','center','VerticalAlignment','middle')
    %     text(count,y_coord+.06,'n.s.','FontSize',25,'HorizontalAlignment','center','VerticalAlignment','middle')
    
%     count = 0;
%     for i_e = o.exp_type_v
%         count = count+1;
%         y_coord = -.5;%0.05;
%         if bias_stat.signif_indiv(i_e) < 0.01
%             text(count,y_coord,'\bf**','FontSize',30,'HorizontalAlignment','center','VerticalAlignment','middle')
%         elseif bias_stat.signif_indiv(i_e) < 0.05
%             text(count,y_coord,'\bf*','FontSize',30,'HorizontalAlignment','center','VerticalAlignment','middle')
%         else
%             text(count,y_coord+.06,'n.s.','FontSize',25,'HorizontalAlignment','center','VerticalAlignment','middle')
%         end
%     end
    
    ax.XAxis.Visible = 'off'
%     plot([1 1 2 2]',[[0 .05 .05  0]+.05]','k','LineWidth',LWS)
%     text(1.5,.15,'n.s.','FontSize',25,'HorizontalAlignment','center','VerticalAlignment','middle')
    
end


% ylabel('Rare $\leftarrow$\textbf{Bias}$\rightarrow$ Freq.','Interpreter','latex')
% set(gca,'FontName','Helvetica','FontSize',20)
% xlabel('valami')
% ylabel('rare \leftarrow \bf bias\rm \rightarrow freq.')

set(gca,'FontName','Calibri Light')
% set(gca,'FontName','Corbel','FontSize',20)
% set(gca,'X')

FS=30;
ax.YAxis.FontSize = FS;
ax.XAxis.FontSize = FS-4;
ax.LineWidth = 2;


ax = gca;               % get the current axis
ax.Clipping = 'off';    % turn clipping off

% set(gca,'position',[.2,.14,.76,.76]);



% fig.PaperUnits = 'centimeters';
% fig.PaperPosition = [0 0 10 5];

if o.plot_save_on
    
    fig_name = ['\figs\bias'];
    for i_e = o.exp_type_v
        fig_name = [fig_name '_e' num2str(i_e)];
    end
    fig_name = [fig_name,'_T',num2str(o.total)]
    curr_dir = pwd;
    % print([curr_dir fig_name],'-dpng','-r600')
    % saveas(fig,[curr_dir fig_name,'.png'])
    
    
    exportgraphics(t,[curr_dir fig_name,'.png'],'Resolution',300)
%     ax.XAxis.Visible = 'off'
    
    
end
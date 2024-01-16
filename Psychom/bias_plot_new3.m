function bias_plot_new3(o, bias_stat)

t = tiledlayout(1,1,'Padding','tight');
t.Units = 'centimeters';
% t.OuterPosition = [1 2 11 10];
if o.plot_num == 1
    t.InnerPosition = [3 2 6 10];
elseif o.plot_num == 2
    t.InnerPosition = [3 2 6 10];
elseif o.plot_num == 10
    t.OuterPosition = [1 2 14 5.5];
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

    % same periods next to each other
    group_width = .8;
    shift_v = linspace(-1,1,length(o.exp_type_v)*2+1)*group_width/2;
    shift_v = shift_v(2:2:end);
    bar_half_width = group_width/length(o.exp_type_v)/2;
    gap_ = .3;
    for i_e = o.exp_type_v
        count=count+1;

        fill([1+shift_v(count)-bar_half_width  1+shift_v(count)+bar_half_width]*[1 1 0 0;0 0 1 1],[0 1 1 0]*bias_stat.mean(i_e,2),o.color(i_e,:),'EdgeColor','none','EdgeColor',o.color(i_e,:)*.5,'LineWidth',2)
        
        hold on
        if ~isnan(bias_stat.sem(i_e,2))
            errorbar(count-bar_half_width/2,bias_stat.mean(i_e,2),bias_stat.sem(i_e,2),'.','Color','k','LineWidth',4)
        end
        
        fill([(1+group_width+gap_)+shift_v(count)-bar_half_width  (1+group_width+gap_)+shift_v(count)+bar_half_width]*[1 1 0 0;0 0 1 1],[0 1 1 0]*bias_stat.mean(i_e,3),o.color(i_e,:),'EdgeColor',o.color(i_e,:)*.5,'LineWidth',2)
        hold on
        if ~isnan(bias_stat.sem(i_e,2))
            errorbar(count-bar_half_width/2,bias_stat.mean(i_e,2),bias_stat.sem(i_e,2),'.','Color','k','LineWidth',4)
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




max_Y = Y_maxpos;


if o.pred


    if o.total | o.group_exp
        axis([.5 count+.5 min_Y max_Y])
        xticks([1:count])



        if o.plot_num == 10 | (o.plot_num == 1 & o.pred_case==1)

            ff=fill([.5*[1 1 1] [1 1 1]*(count+.5)],[-1 0 1 1 0 -1]*max(abs(min_Y),abs(max_Y))/6,[1 1 1]*.5,'EdgeColor','none')
            ff.FaceVertexAlphaData = [0 1 0 0 1 0]';
            ff.FaceAlpha = 'interp';

        else

            plot([.5 count+.5],[0 0],'Color',[1 1 1]*.15,'LineWidth',3)

        end

    else
        axis([.5 (1+group_width+gap_)+.5  min_Y max_Y])
        %     xticks([1 2])
        xticks([])

        %plot([.5 (1+group_width+gap_)+.5],[0 0],'k','LineWidth',3)
        [.5 (1+group_width+gap_)+.5 (1+group_width+gap_)+.5 .5]
        [-1 1 1 -1]*max(abs(min_Y),abs(max_Y))/10
        %     plot([.5 (1+group_width+gap_)+.5],[0 0],'k','LineWidth',3)

        if o.pred

            switch o.pred_case
                case {3,4}
                    ff=fill([.5*[1 1 1] [1 1 1]*(1+group_width+gap_)+.5],[-1 0 1 1 0 -1]*max(abs(min_Y),abs(max_Y))/6,[1 1 1]*.5,'EdgeColor','none')
                    ff.FaceVertexAlphaData = [0 1 0 0 1 0]';
                    %     ff=fill([.5*[1 1] [1 1]*(1+group_width+gap_)+.5],[0 1 1 0]*max(abs(min_Y),abs(max_Y))/8,[1 1 1]*.5,'EdgeColor','none')
                    %     ff.FaceVertexAlphaData = [1 0 0 1]';
                    ff.FaceAlpha = 'interp';
                    'x'
                case {1,2}
                    plot([.5 (1+group_width+gap_)+.5],[0 0],'Color',[1 1 1]*.15,'LineWidth',3)
            end
        end

    end
    % xticklabels(exp_name(o.exp_type_v))
    xticklabels([])
    % yticks([-1:.1:1])
    % yticks([-1.5:.3:1.5])
    if o.pred
        yticks([])
    else
        yticks([-1.2:.4:1.2])
    end

    box off

end



ax = gca;
if o.plot_num == 1


   axis([.5 4.5 -.3 .7])
    
   bias_stat.signif_indiv
   count = 0;
    
   o.color_ND = [0.8320  0.1563  0.0195];



   xticklabels(exp_name(o.exp_type_v))
   ax.XAxis.Visible = 'off';
   
elseif o.plot_num == 2
   
    axis([.5 3.5 -.25 .15])
    xticklabels(exp_name(o.exp_type_v))
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

    switch o.pred_case
        case 1

            for k = 1:2
                plot([1 1 1 1]*(1+(k-1)*(gap_+group_width))'+[-1 -1 2 2]*group_width/6,[[0 .03 .03  -.14]+.84]','k','LineWidth',LWS)
                text((1+(k-1)*(gap_+group_width))+group_width/6/2,.89,'\bf*','FontSize',35,'HorizontalAlignment','center','VerticalAlignment','middle','Color',[0 0 0])
                plot([1 1 1 1]*(1+(k-1)*(gap_+group_width))'+[-2 -2 0 0]*group_width/6,[[0 .03 .03  0]+.70]','k','LineWidth',LWS)
                text((1+(k-1)*(gap_+group_width))-group_width/6,.8,'\bf n.s.','FontSize',25,'HorizontalAlignment','center','VerticalAlignment','middle','Color',[0 0 0])
            end
            yticks([0])
            yticklabels({'\bf 0'})

        case 2
            ax.XAxisLocation = 'origin';
            for k = 1:2
                plot([1 1 1 1]*((k-1)*(gap_+group_width)+1)'+[-1 -1 2 2]*group_width/6,[[0 .03 .03  -.14]+.18]','k','LineWidth',LWS)
                text((1+(k-1)*(gap_+group_width))+group_width/6/2,.23,'\bf*','FontSize',35,'HorizontalAlignment','center','VerticalAlignment','middle','Color',[0 0 0])
                plot([1 1 1 1]*(1+(k-1)*(gap_+group_width))'+[-2 -2 0 0]*group_width/6,[[0 .03 .03  0]+.04]','k','LineWidth',LWS)
                text((1+(k-1)*(gap_+group_width))-group_width/6,.14,'\bf n.s.','FontSize',25,'HorizontalAlignment','center','VerticalAlignment','middle','Color',[0 0 0])
            end

            yticks([0])

            yticklabels({'\bf 0'})

        case 3
            ax.XAxis.Visible = 'off';

            for k = 1:2
                plot([1 1 1 1]*(1+(k-1)*(gap_+group_width))'+[-2 -2 -.2 -.2]*group_width/6,[[0 .05 .05  0]+.55]','k','LineWidth',LWS)
                plot([1 1 1 1]*(1+(k-1)*(gap_+group_width))'+[.2 .2   2   2]*group_width/6,[[0 .05 .05  0]+.55]','k','LineWidth',LWS)
                text((1+(k-1)*(gap_+group_width))-group_width/6*1.1,.58+.05,'\bf*','FontSize',35,'HorizontalAlignment','center','VerticalAlignment','middle','Color',[0 0 0])
                text((1+(k-1)*(gap_+group_width))+group_width/6*1.1,.58+.05,'\bf*','FontSize',35,'HorizontalAlignment','center','VerticalAlignment','middle','Color',[0 0 0])
            end

        case 4

            ax.XAxis.Visible = 'off';
            for k = 1:2
                shift = 0;
                plot([1 1 1 1]*(1+(k-1)*(gap_+group_width))'+[-2 -2 -.2 -.2]*group_width/6,[[0 .05 .05  0]+.55]'+shift,'k','LineWidth',LWS)
                plot([1 1 1 1]*(1+(k-1)*(gap_+group_width))'+[.2 .2   2   2]*group_width/6,[[0 .05 .05  0]+.55]'+shift,'k','LineWidth',LWS)
                text((1+(k-1)*(gap_+group_width))-group_width/6*1.1,.58+.05+shift,'\bf*','FontSize',35,'HorizontalAlignment','center','VerticalAlignment','middle','Color',[0 0 0])
                text((1+(k-1)*(gap_+group_width))+group_width/6*1.1,.58+.05+shift,'\bf*','FontSize',35,'HorizontalAlignment','center','VerticalAlignment','middle','Color',[0 0 0])
            end

    end


end




set(gca,'FontName','Calibri Light')
% set(gca,'FontName','Corbel','FontSize',20)
% set(gca,'X')

FS=30;
ax.YAxis.FontSize = FS;
ax.XAxis.FontSize = FS-4;
ax.LineWidth = 3;


ax = gca;               % get the current axis
ax.Clipping = 'off';    % turn clipping off


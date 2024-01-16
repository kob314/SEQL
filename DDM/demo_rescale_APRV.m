%%
clear all
close all

set(0,'DefaultFigureWindowStyle','docked')

%% load in data

include_General

o.exp_type_v = [1 6 7 4 3 2 14 17];
o.color = [0.4 0.55 0.25; .45 .25 .5;.9 .35 .4; .2 .3 .6; 0 0 0; .86 .45 .1; 0.93 0.74 0.35];
o.color([14 17],:) = [.6 .35 .2; [1 1 1]*.3];
true_AP = [.65 .5 .65 .65 0 .65 .5 0 0 0 0 0 0 .65 0 0 .65];

%%
load("DATA_EXP4.mat")
subj_AP = cell2mat(arrayfun(@(ie) arrayfun(@(is) DATA_EXP(ie).subj(is).subj_AP, DATA_EXP(ie).incl_subj), o.exp_type_v,'UniformOutput',false))';

t=tiledlayout(2,1,'Units','pixels','InnerPosition',[10 10 1000 1000]);
nexttile(t)
hold on

tt    = linspace(0,2*pi,100);
sin_t = sin(tt);
cos_t = cos(tt);

o.exp_type_v = [7 1 4 6 3 2];
% o.exp_type_v = [1 4 6 3];
for j = 1:2

    % res0 = load('data_numpyro_static_new_subjAP_ysc2.mat');
    % % res0 = load('data_numpyro_static_new_subjAP_ysc2_ym3.mat');
    % % res0 = load('numpyro_static_new (4).mat');
    % 
    % switch j
    %     case 1
    %         res1 = load('data_numpyro_static_new_subjAP_ysc1.mat');
    %     case 2
    %         res1 = load('data_numpyro_static_new_subjAP_ysc15.mat')
    %     case 3
    %         res1 = load('data_numpyro_static_new_subjAP_ysc25.mat');
    %     case 4
    %         res1 = load('data_numpyro_static_new_subjAP_ysc3.mat');
    % end


    res0 = load('data_numpyro_static_new_subjAP_ysc2_ym3_simpleAP.mat');
    switch j
        case 1
            res1 = load('data_numpyro_static_new_subjAP_ysc1_ym3_simpleAP.mat');
        case 2
            res1 = load('data_numpyro_static_new_subjAP_ysc3_ym3_simpleAP.mat');
    end

    for i = 1:2
        switch i
            case 1
                res = res0;
            case 2
                res = res1;
        end

        op_pair_effect = 0;
        aux_q = res.subj.q(:,1);
        q_M = accumarray(res.subj_struct(:,1)+1, res.subj.q(:,1) - op_pair_effect,[],@mean)';
        q_L = accumarray(res.subj_struct(:,1)+1, res.subj.q(:,1),[],@length)';
        q_S = accumarray(res.subj_struct(:,1)+1, res.subj.q(:,1)- op_pair_effect,[],@std)'./sqrt(q_L);

        aux_nu = subj_AP - res.subj.q(:,1) +.5;
        nu_M = accumarray(res.subj_struct(:,1)+1, aux_nu,[],@mean)';
        nu_L = accumarray(res.subj_struct(:,1)+1, aux_nu,[],@length)';
        nu_S = accumarray(res.subj_struct(:,1)+1, aux_nu,[],@std)'./sqrt(q_L);

        q_(i).M = q_M; q_(i).S = q_S;
        nu_(i).M = nu_M; nu_(i).S = nu_S;
    end

    plot(q_(1).M(o.exp_type_v),q_(2).M(o.exp_type_v),'g','LineWidth',2)
     % plot(q_(1).M(o.exp_type_v),(q_(1).M(o.exp_type_v)-.5)*1.5+.5,'g','LineWidth',2)
    scatter(q_(1).M(o.exp_type_v),q_(2).M(o.exp_type_v),180,o.color(o.exp_type_v,:),'filled')
    for ie = o.exp_type_v
        fill(q_(1).M(ie)+q_(1).S(ie).*sin_t',q_(2).M(ie)+q_(2).S(ie).*cos_t',o.color(ie,:),'FaceAlpha',.2,'EdgeColor','none')%,o.color(ie,:),'EdgeAlpha',.2)
    end
    plot([.4 .7],[.4 .7],'k')


end
xlabel('AP')
xlabel('RV')
xticks(.1:.1:.9)
yticks(.1:.1:.9)
% d=.11;
d=.09;
d=.25;
axis([-d d -d d]+.535)
axis square
set(gca,'LineWidth',2,'FontSize',26,'FontName','Calibri Light')
box off


%%
return

nexttile(t)
hold on
scatter(q_(1).M(o.exp_type_v),q_(2).M(o.exp_type_v),180,o.color(o.exp_type_v,:),'filled')
for ie = o.exp_type_v
    ie
    fill(q_(1).M(ie)+q_(1).S(ie).*sin_t',q_(2).M(ie)+q_(2).S(ie).*cos_t',o.color(ie,:),'FaceAlpha',.5,'EdgeColor',o.color(ie,:),'EdgeAlpha',.9)
    % plot([q_(1).M(ie)-q_(1).S(ie); q_(1).M(ie)+q_(1).S(ie)],[q_(2).M(ie)].*[1 1]','color',o.color(ie,:),'LineWidth',2)
    % plot([q_(1).M(ie)].*[1 1]',[q_(2).M(ie)-q_(2).S(ie); q_(2).M(ie)+q_(2).S(ie)],'color',o.color(ie,:),'LineWidth',2)
    % plot([q_(1).M(ie)-q_(1).S(ie); q_(1).M(ie)+q_(1).S(ie)],[q_(2).M(ie)].*[1 1]','color',o.color(ie,:),'LineWidth',2)
end
plot([.4 .7],[.4 .7],'k')

xticks(.1:.2:.9)
yticks(.1:.2:.9)
d=.15;
axis([-d d -d d]+.55)
axis square
set(gca,'LineWidth',2,'FontSize',26,'FontName','Calibri Light')
box off

nexttile(t)
hold on
plot([0 1; .5 .5]',[.5 .5; 0 1]','color',[1 1 1]*0,'LineWidth',1)
plot([1 0]'+[0 .15],[0 1]'.*[1 1],'color',[1 1 1]*.65,'LineWidth',3)



tt    = linspace(0,2*pi,100);
sin_t = sin(tt);
cos_t = cos(tt);

jit_v = [-1 0 -1 -1.3 0 1.3 0 0 0 0 0 0 0 1.3 0 0 1]*.006;
i=2
for ie = o.exp_type_v
    % fill(q_M(ie)+q_S(ie).*sin_t',nu_M(ie)+q_S(ie).*cos_t',o.color(ie,:),'FaceAlpha',.5,'EdgeColor',o.color(ie,:),'EdgeAlpha',.9)
    % plot(q_M(ie)+jit_v(ie)+q_S(ie).*[sin(pi/4) -cos(pi/4)],nu_M(ie)+jit_v(ie)+q_S(ie).*[-sin(pi/4) cos(pi/4)],'LineWidth',5,'Color',o.color(ie,:))
    rotate_errorbar([q_(i).M(ie) nu_(i).M(ie)]+jit_v(ie),q_(i).S(ie),-pi/4,o.color(ie,:),4,.01)

    % scatter(aux_q(res.subj_struct(:,1)+1==ie),aux_nu(res.subj_struct(:,1)+1==ie),10,o.color(ie,:),"filled")
end
scatter(q_(i).M(o.exp_type_v)+jit_v(o.exp_type_v),nu_(i).M(o.exp_type_v)+jit_v(o.exp_type_v),180,o.color(o.exp_type_v,:),'filled')%,'MarkerEdgeColor',[1 1 1]*.3,'LineWidth',2)
xticks([.1:.2:.9])
yticks([.1:.2:.9])
% xticks([.5:.15:.9])
% yticks([.5:.15:.9])
d=.2

d=1
axis([-d d -d d]+.5)
% axis([.42 .68 .42 .68])
axis square

set(gca,'LineWidth',2,'FontSize',26,'FontName','Calibri Light')
box off



% %% old
% 
% 
% %%
% clear all
% close all
% 
% set(0,'DefaultFigureWindowStyle','docked')
% 
% %% load in data
% 
% include_General
% % data_ = data_reader;
% 
% res0 = load('data_numpyro_static_new_subjAP_ysc2.mat');
% res1 = load('data_numpyro_static_new_subjAP_ysc25.mat');
% % res1 = load('data_numpyro_static_new_subjAP_ysc3.mat');
% % res1 = load('data_numpyro_static_new_subjAP_ysc2_ym2.mat');
% 
% 
% o.exp_type_v = [1 6 7 4 3 2 14 17];
% o.color = [0.4 0.55 0.25; .45 .25 .5;.9 .35 .4; .2 .3 .6; 0 0 0; .86 .45 .1; 0.93 0.74 0.35];
% o.color([14 17],:) = [.6 .35 .2; [1 1 1]*.3];
% true_AP = [.65 .5 .65 .65 0 .65 .5 0 0 0 0 0 0 .65 0 0 .65];
% 
% %%
% load("DATA_EXP4.mat")
% subj_AP = cell2mat(arrayfun(@(ie) arrayfun(@(is) DATA_EXP(ie).subj(is).subj_AP, DATA_EXP(ie).incl_subj), o.exp_type_v,'UniformOutput',false))';
% 
% for i = 1:2
%     switch i
%         case 1
%             res = res0;
%         case 2
%             res = res1;
%     end
% 
%     op_pair_effect = 0;
%     aux_q = res.subj.q(:,1);
%     q_M = accumarray(res.subj_struct(:,1)+1, res.subj.q(:,1) - op_pair_effect,[],@mean)';
%     q_L = accumarray(res.subj_struct(:,1)+1, res.subj.q(:,1),[],@length)';
%     q_S = accumarray(res.subj_struct(:,1)+1, res.subj.q(:,1)- op_pair_effect,[],@std)'./sqrt(q_L);
% 
%     aux_nu = subj_AP - res.subj.q(:,1) +.5;
%     nu_M = accumarray(res.subj_struct(:,1)+1, aux_nu,[],@mean)';
%     nu_L = accumarray(res.subj_struct(:,1)+1, aux_nu,[],@length)';
%     nu_S = accumarray(res.subj_struct(:,1)+1, aux_nu,[],@std)'./sqrt(q_L);
% 
%     q_(i).M = q_M; q_(i).S = q_S;
%     nu_(i).M = nu_M; nu_(i).S = nu_S;
% end
% 
% t=tiledlayout(2,1,'Units','pixels','InnerPosition',[100 100 400 800])
% 
% tt    = linspace(0,2*pi,100);
% sin_t = sin(tt);
% cos_t = cos(tt);
% 
% nexttile(t)
% hold on
% scatter(q_(1).M(o.exp_type_v),q_(2).M(o.exp_type_v),180,o.color(o.exp_type_v,:),'filled')
% for ie = o.exp_type_v(1:8)
%     ie
%     fill(q_(1).M(ie)+q_(1).S(ie).*sin_t',q_(2).M(ie)+q_(2).S(ie).*cos_t',o.color(ie,:),'FaceAlpha',.5,'EdgeColor',o.color(ie,:),'EdgeAlpha',.9)
%     % plot([q_(1).M(ie)-q_(1).S(ie); q_(1).M(ie)+q_(1).S(ie)],[q_(2).M(ie)].*[1 1]','color',o.color(ie,:),'LineWidth',2)
%     % plot([q_(1).M(ie)].*[1 1]',[q_(2).M(ie)-q_(2).S(ie); q_(2).M(ie)+q_(2).S(ie)],'color',o.color(ie,:),'LineWidth',2)
%     % plot([q_(1).M(ie)-q_(1).S(ie); q_(1).M(ie)+q_(1).S(ie)],[q_(2).M(ie)].*[1 1]','color',o.color(ie,:),'LineWidth',2)
% end
% plot([.4 .7],[.4 .7],'k')
% 
% xticks(.1:.2:.9)
% yticks(.1:.2:.9)
% d=.15;
% axis([-d d -d d]+.55)
% axis square
% set(gca,'LineWidth',2,'FontSize',26,'FontName','Calibri Light')
% box off
% 
% nexttile(t)
% hold on
% plot([0 1; .5 .5]',[.5 .5; 0 1]','color',[1 1 1]*0,'LineWidth',1)
% plot([1 0]'+[0 .15],[0 1]'.*[1 1],'color',[1 1 1]*.65,'LineWidth',3)
% 
% 
% 
% tt    = linspace(0,2*pi,100);
% sin_t = sin(tt);
% cos_t = cos(tt);
% 
% jit_v = [-1 0 -1 -1.3 0 1.3 0 0 0 0 0 0 0 1.3 0 0 1]*.006;
% i=2
% for ie = o.exp_type_v
%     % fill(q_M(ie)+q_S(ie).*sin_t',nu_M(ie)+q_S(ie).*cos_t',o.color(ie,:),'FaceAlpha',.5,'EdgeColor',o.color(ie,:),'EdgeAlpha',.9)
%     % plot(q_M(ie)+jit_v(ie)+q_S(ie).*[sin(pi/4) -cos(pi/4)],nu_M(ie)+jit_v(ie)+q_S(ie).*[-sin(pi/4) cos(pi/4)],'LineWidth',5,'Color',o.color(ie,:))
%     rotate_errorbar([q_(i).M(ie) nu_(i).M(ie)]+jit_v(ie),q_(i).S(ie),-pi/4,o.color(ie,:),4,.01)
% 
%     % scatter(aux_q(res.subj_struct(:,1)+1==ie),aux_nu(res.subj_struct(:,1)+1==ie),10,o.color(ie,:),"filled")
% end
% scatter(q_(i).M(o.exp_type_v)+jit_v(o.exp_type_v),nu_(i).M(o.exp_type_v)+jit_v(o.exp_type_v),180,o.color(o.exp_type_v,:),'filled')%,'MarkerEdgeColor',[1 1 1]*.3,'LineWidth',2)
% xticks([.1:.2:.9])
% yticks([.1:.2:.9])
% % xticks([.5:.15:.9])
% % yticks([.5:.15:.9])
% d=.2
% axis([-d d -d d]+.5)
% % axis([.42 .68 .42 .68])
% axis square
% 
% set(gca,'LineWidth',2,'FontSize',26,'FontName','Calibri Light')
% box off

clear all
close all
set(0,'DefaultFigureWindowStyle','docked')

% load...
% res = load("pymc_static_fit.mat")
% res = load("pymc_static_fit_fixTau.mat")
% res = load("numpyro_static_new.mat")
res = load("numpyro_static_new (4).mat");
% res = load("numpyro_static_new_subjAP.mat");



o.exp_type_v = [1 6 7 4 3 2 14 17];
o.color = [0.4 0.55 0.25; .45 .25 .5;.9 .35 .4; .2 .3 .6; 0 0 0; .86 .45 .1; 0.93 0.74 0.35];
o.color([14 17],:) = [.6 .35 .2; [1 1 1]*.3];
% true_AP = [.65 .5 .65 .65 0 .65 .5];
true_AP = [.65 .5 .65 .65 0 .65 .5 0 0 0 0 0 0 .65 0 0 .65];

%%
load("DATA_EXP4.mat")
load('obj_correction_1_2_3_4_6_7_14_17.mat')

ob_pair_m = [];
for ie = o.exp_type_v
    ob_pair = arrayfun(@(is) [DATA_EXP(ie).subj(is).obj_pair], DATA_EXP(ie).incl_subj,'UniformOutput',false);
    ob_pair = cell2mat(ob_pair');
    ob_pair_m = [ob_pair_m; ob_pair];
end
op_pair_effect = obj_bias(ob_pair_m(:,2)) - obj_bias(ob_pair_m(:,1)); 
op_pair_effect = op_pair_effect*0;


subj_AP = cell2mat(arrayfun(@(ie) arrayfun(@(is) DATA_EXP(ie).subj(is).subj_AP, DATA_EXP(ie).incl_subj), o.exp_type_v,'UniformOutput',false))';

aux_q = res.subj.q(:,1);
q_M = accumarray(res.subj_struct(:,1)+1, res.subj.q(:,1) - op_pair_effect,[],@mean)';
q_L = accumarray(res.subj_struct(:,1)+1, res.subj.q(:,1),[],@length)';
q_S = accumarray(res.subj_struct(:,1)+1, res.subj.q(:,1)- op_pair_effect,[],@std)'./sqrt(q_L);

aux_nu = subj_AP - res.subj.q(:,1) +.5;
nu_M = accumarray(res.subj_struct(:,1)+1, aux_nu,[],@mean)';
nu_L = accumarray(res.subj_struct(:,1)+1, aux_nu,[],@length)';
nu_S = accumarray(res.subj_struct(:,1)+1, aux_nu,[],@std)'./sqrt(q_L);

% q_M_obj_bias = accumarray(res.subj_struct(:,1)+1, - op_pair_effect,[],@mean)';
% q_M(o.exp_type_v) = 1./(1+exp( -1.5*res.exp.q_exp(:,1) ))';
% q_S(o.exp_type_v) = [1./(1+exp( -1.5*( res.exp.q_exp(:,1)+res.exp.q_exp(:,2) ) ))-1./(1+exp(-1.5*res.exp.q_exp(:,1) ))]';

% % q_M = q_M + q_M_obj_bias

%%

% nu_M = true_AP-q_M+.5;

% o.exp_type_v = [1 6 7 14 17];
t=tiledlayout(1,1,'Units','pixels','InnerPosition',[100 100 500 500])

nexttile(t)
hold on
plot([0 1; .5 .5]',[.5 .5; 0 1]','color',[1 1 1]*0,'LineWidth',1)
plot([1 0]'+[0 .15],[0 1]'.*[1 1],'color',[1 1 1]*.65,'LineWidth',3)


tt    = linspace(0,2*pi,100);
sin_t = sin(tt);
cos_t = cos(tt);

jit_v = [-1 0 -1 -1.3 0 1.3 0 0 0 0 0 0 0 1.3 0 0 1]*.006;
for ie = o.exp_type_v
    % fill(q_M(ie)+q_S(ie).*sin_t',nu_M(ie)+q_S(ie).*cos_t',o.color(ie,:),'FaceAlpha',.5,'EdgeColor',o.color(ie,:),'EdgeAlpha',.9)
    % plot(q_M(ie)+jit_v(ie)+q_S(ie).*[sin(pi/4) -cos(pi/4)],nu_M(ie)+jit_v(ie)+q_S(ie).*[-sin(pi/4) cos(pi/4)],'LineWidth',5,'Color',o.color(ie,:))
    rotate_errorbar([q_M(ie) nu_M(ie)]+jit_v(ie),q_S(ie),-pi/4,o.color(ie,:),4,.01)

    % scatter(aux_q(res.subj_struct(:,1)+1==ie),aux_nu(res.subj_struct(:,1)+1==ie),10,o.color(ie,:),"filled")
end
scatter(q_M(o.exp_type_v)+jit_v(o.exp_type_v),nu_M(o.exp_type_v)+jit_v(o.exp_type_v),180,o.color(o.exp_type_v,:),'filled')%,'MarkerEdgeColor',[1 1 1]*.3,'LineWidth',2)
xticks([.1:.2:.9])
yticks([.1:.2:.9])
% xticks([.5:.15:.9])
% yticks([.5:.15:.9])
d=.2
axis([-d d -d d]+.5)
% axis([.42 .68 .42 .68])
axis square

set(gca,'LineWidth',2,'FontSize',26,'FontName','Calibri Light')
box off

exportgraphics(t,'figures\demo_statBayes_params.emf','ContentType','vector')
% exportgraphics(t,'figures\demo_statBayes_params_poster.emf','ContentType','vector')
rt
%%
close

t=tiledlayout(2,1,'Units','pixels','InnerPosition',[100 100 300 500])

nexttile(t)
hold on

d=.1
count = 0;
for ie = o.exp_type_v
    count = count+1;
    fill(count+[d d 1-d 1-d],(q_M(ie)-.5)*[0 1 1 0]+.5,o.color(ie,:),'FaceAlpha',1,'EdgeColor',o.color(ie,:)*.7,'linewidth',3)
end

plot([1 count+1],[.5 .5],'k','LineWidth',3)
xlim([1 count+1])

axis([1 count+1 .3 .7])
yticks([.1:.2:.9])
xticks([])
ax = gca;
ax.XAxis.Visible = 'off';

set(gca,'LineWidth',2,'FontSize',26,'FontName','Calibri Light')
box off



nexttile(t)
hold on

count = 0;
for ie = o.exp_type_v
    count = count+1;
    fill(count+[d d 1-d 1-d],(nu_M(ie)-.5)*[0 1 1 0]+.5,o.color(ie,:),'FaceAlpha',1,'EdgeColor',o.color(ie,:)*.7,'linewidth',3)
end

plot([1 count+1],[.5 .5],'k','LineWidth',3)

axis([1 count+1 .3 .7])
yticks([.1:.2:.9])
xticks([])
ax = gca;
ax.XAxis.Visible = 'off';

% xticks(1:6); xticklabels({'E1'})

set(gca,'LineWidth',2,'FontSize',26,'FontName','Calibri Light')
box off

exportgraphics(t,'figures\demo_statBayes_param_biases.emf','ContentType','vector')

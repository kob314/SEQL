clear all
close all
set(0,'DefaultFigureWindowStyle','docked')

%% grids
o.N_grid.x = 100;
o.N_grid.y = 80;
o.N_grid.nu = 41;
o.N_grid.q = 41;


%%
test_AP = .65;


par_info = struct( ...
    'q' , 1, ...
    'a' , 1, ...
    'b' , 1, ...
    'c' , 1, ...
    'beta'    , 1, ...
    'lambda'  , 1, ...
    'kappa'   , 1, ...
    'tau'     , 1);

par = struct( ...
    'q' , .5, ...
    'a' , 1, ...
    'b' , 1, ...
    'c' , 1, ...
    'beta'    , 1, ...
    'lambda'  , 1, ...
    'kappa'   , 1, ...
    'tau'     , 1);


par_v = parameter_wrap(par,par_info,'s2v');
par_s = parameter_wrap(par_v,par_info,'v2s');

data_in.y_min_max = [.025 1];
par_s.nu = .5 + (test_AP - par_s.q);
distr = pre_calc_data_generator_static( par_s , data_in ,o);

col_ = lines(5);
cols_ = [0 176 240; 213 40 5]/255;
LW = 5;


close
t = tiledlayout(4,2,'Units','centimeters','InnerPosition',[2 3 20 26]);


for j = 1:2
    switch j
        case 1
            par_s.q = .65;
        case 2
            par_s.q = .5;
    end
    par_s.nu = .5 + (test_AP - par_s.q);
    distr = pre_calc_data_generator_static( par_s , data_in ,o);
    
    for i=1:4
        nexttile(2*(i-1)+j)
        hold on
        switch i
            case 1
                plot([0 200 200 500],[0.5 0.5 0.65 0.65],'LineWidth',LW,'Color',cols_(j,:))
                axis([0 500 .35 .8])
                xlabel('time'); xticks([0 200 500])
                if j==1; ylabel('$AP$','Interpreter','latex'); else; ylabel('$ND$','Interpreter','latex'); end
            case 2
                fill([.2 .2 .8 .8]  ,[0 1 1 0]*exp(distr.lP_z_G_q(:,:,1)),cols_(1,:)*.1+.9*[1 1 1],'EdgeColor',cols_(1,:),'LineWidth',LW,'LineStyle',':')
                fill([.2 .2 .8 .8]+1,[0 1 1 0]*exp(distr.lP_z_G_q(:,:,2)),cols_(1,:)*.1+.9*[1 1 1],'EdgeColor',cols_(1,:),'LineWidth',LW,'LineStyle','-')
                axis([0 2 0 1])
                xticks([0 1]+.5); xticklabels({'rare','frequent'})
                xlabel('$z$','Interpreter','latex');
                ylabel('$P(z)$','Interpreter','latex');
            case 3
                ph(1) = plot(distr.grid.x,exp(distr.lP_x_G_z_nu(:,:,1)),':','LineWidth',LW,'Color',cols_(2,:));
                ph(2) = plot(distr.grid.x,exp(distr.lP_x_G_z_nu(:,:,2)),'-','LineWidth',LW,'Color',cols_(2,:));
                xlabel('$x$','Interpreter','latex'); xlim([-4 4])
                ylabel('$P(x|z)$','Interpreter','latex'); yticks([])
            case 4
                ph(1) = plot(distr.grid.x,exp(distr.lP_x_G_z_nu(:,:,1)+distr.lP_z_G_q(:,:,1)),':','Color',[1 1 1]*.5,'LineWidth',LW);
                ph(2) = plot(distr.grid.x,exp(distr.lP_x_G_z_nu(:,:,2)+distr.lP_z_G_q(:,:,2)),'-','Color',[1 1 1]*.5,'LineWidth',LW);
                if j==2
                    plot(distr.grid.x,exp(distr.lP_x_G_z_nu(:,:,1)+distr.lP_z_G_q(:,:,1)) + exp(distr.lP_x_G_z_nu(:,:,2)+distr.lP_z_G_q(:,:,2)),'color',cols_(1,:),'LineWidth',LW*2,'LineStyle','--');
                end
                ph(3) = plot(distr.grid.x,exp(distr.lP_x_G_z_nu(:,:,1)+distr.lP_z_G_q(:,:,1)) + exp(distr.lP_x_G_z_nu(:,:,2)+distr.lP_z_G_q(:,:,2)),'color',cols_(j,:),'LineWidth',LW);
                xlabel('$x$','Interpreter','latex'); xlim([-4 4])
                ylabel('$P(x)$','Interpreter','latex'); yticks([])
        end
        yticks([])
        set(gca,'LineWidth',2,'FontSize',26,'FontName','Calibri Light')
        box off
    end
end

% exportgraphics(t,'figures\two_interp_demo.emf')
% print('figures\two_interp_demo.svg','-e')


% legend(ph,{'$P(x,z=\mathrm{rare})$','$P(x,z=\mathrm{freq.})$','$P(x)$'},'Interpreter','latex'); legend boxoff

% plot(distr.grid.x,exp(distr.lP_x_G_z_unbiased(:,:,1)+distr.lP_z_G_q(:,:,1)),'--k','LineWidth',LW)
% hold on
% plot(distr.grid.x,exp(distr.lP_x_G_z_unbiased(:,:,2)+distr.lP_z_G_q(:,:,2)), '-k','LineWidth',LW)
% plot(distr.grid.x,exp(distr.lP_x_G_z_unbiased(:,:,1)+distr.lP_z_G_q(:,:,1)) + exp(distr.lP_x_G_z_unbiased(:,:,2)+distr.lP_z_G_q(:,:,2)),'color',lines(1),'LineWidth',LW)




% xticks([-4 0 4])

% xlabel('$x$','Interpreter','latex'); xlim([-4 4])
% ylabel('$P(x)$','Interpreter','latex'); yticks([])
% set(gca,'LineWidth',2,'FontSize',26,'FontName','Calibri Light')
% box off


% plot(distr.grid.x,exp(distr.lP_x_G_z_unbiased(:,:,1)+log(1-test_AP))+exp(distr.lP_x_G_z_unbiased(:,:,2)+log(test_AP)),'Color',col_(1,:),'LineWidth',3)
    
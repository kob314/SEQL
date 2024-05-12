clear all
close all

set(0,'DefaultFigureWindowStyle','docked')

addpath('customcolormap\')

%% choose plot
%%% 5 DIFFERENT PLOTS

o.plot=4 %1: LLH surface; 2: predictive distributions; 3: psyhcometrics  %% 4 BIAS surf %% Bias 2 % bias without contours

%% settings
% colors & sizes
color_AP = [0       0.6875  0.9375];
color_ND = [0.8320  0.1563  0.0195];

FS = 40;
MS = 1200;
LW = 14;
LW_UB = 10;

%% model parameters
par_info = struct( ...
    'a' , 1, ...
    'b' , 1, ...
    'c' , 1, ...
    'beta'    , 1, ...
    'lambda'  , 1, ...
    'kappa'   , 1, ...
    'tau'     , 1, ...
    'dyn_conc_nu', 1,...
    'dyn_conc_q',  1,...
    'mix_nu', 1,...
    'mix_q',  1);


par = struct( ...
    'a' , .55, ...
    'b' , .4, ...
    'c' ,  3, ...
    'beta'    , 20, ...
    'lambda'  , 1, ...
    'kappa'   , 1, ...
    'tau'     , 1, ...
    'dyn_conc_nu', 1,...
    'dyn_conc_q',  1,...
    'mix_nu', 1,...
    'mix_q',  1);


par_v = parameter_wrap(par,par_info,'s2v');
par_s = parameter_wrap(par_v,par_info,'v2s');

%% grids
distr.N_grid.x = 100;
distr.N_grid.y = 80;
distr.N_grid.nu = 81;
distr.N_grid.q  = 81;

distr.y_min_max = [.045 1];

%% compute the distribution
distr = pre_calc_data_generator_dynamic( par_s , distr);

idx_nu = 41;
idx_q  = 54;

disp(['nu: ',num2str(distr.grid.nu(idx_nu)),'  q: ',num2str(distr.grid.q(idx_q))])

aux=(exp( 10*(sum(log( distr.P_x_G_nu_q(:,:,:,idx_nu,idx_q) ) .* distr.P_x_G_nu_q(:,:,:,idx_nu,idx_q)) - sum(log( distr.P_x_G_nu_q(:,:,:,idx_nu,idx_nu) ) .* distr.P_x_G_nu_q(:,:,:,idx_nu,idx_q)) )) -1 )*100;
disp(['best solution is ', num2str(round(aux,2)),'% better then the unbiased solution'])

avglLH_diff = squeeze(sum(log(distr.P_x_G_nu_q) .* distr.P_x_G_nu_q(:,:,:,idx_nu,idx_q)));
KL_div = squeeze( sum( (log(distr.P_x_G_nu_q(:,:,:,idx_nu,idx_q)) - log(distr.P_x_G_nu_q)).*distr.P_x_G_nu_q(:,:,:,idx_nu,idx_q) ,1) );

%% plots
close

t = tiledlayout(1,1,'TileSpacing','compact');
t.Units = 'pixels';

%%% crop settings
crop         = [.45 .7];
disp(['crop: ',mat2str(crop)])
[~,crop_idx] = min(abs(distr.grid.q(:)-crop));
crop_idx     = crop_idx + [-1,1];
crop_idx = min( crop_idx, length(distr.grid.q)); crop_idx = max( crop_idx, 1);
grid_q  = squeeze(distr.grid.q(crop_idx(1):crop_idx(2)));
grid_nu = squeeze(distr.grid.nu(crop_idx(1):crop_idx(2)));

if o.plot == 1 %% likelihood surface

    t.InnerPosition = [120 120 500 500];
    nexttile(t)

    Z = avglLH_diff(crop_idx(1):crop_idx(2),crop_idx(1):crop_idx(2));

    scale = 3;
    LHcolor = ([1 .8 0]*.7+.3*[1 1 1])*.9;
    mymap = linspace(0,1,1024)'.^scale.*LHcolor(1,:) + (1-linspace(0,1,1024)'.^scale).*ones(1,3);
    colormap(mymap);

    pch = pcolor(grid_q,grid_nu,Z);
    view(0,90)
    % colorbar('box','off')
    shading interp
    hold on


    contour(grid_q,grid_nu,Z,min(Z(:))+(max(Z(:))-min(Z(:)))*[0:.2:.9 .995],'LineColor',LHcolor(1,:)*.8+[1 1 1]*0,'LineWidth',3)%,'ShowText','on')
    view(0,90)


    scatter3([distr.grid.q(idx_q)],[distr.grid.nu(idx_nu)],[.1],MS,color_AP,'filled','LineWidth',2,'MarkerEdgeColor',color_AP*.6)
    scatter3([distr.grid.nu(idx_nu)],[distr.grid.q(idx_q)],[.1],MS,color_ND,'filled','LineWidth',2,'MarkerEdgeColor',color_ND*.6)
    scatter3([.5],[.5],[.1],MS,'k','LineWidth',10,'Marker','+')


    axis([crop crop])
    xticks([.35 .5 .65])
    yticks([.35 .5 .65])


    % title('parameter likelihood')
    box off
    set(gca,'FontSize',FS,'FontName','Calibri Light','layer','top','LineWidth',2)

    exportgraphics(t,'figures\static_LH_surface.emf','ContentType','vector')
%     print('figures\static_LH_surface.svg')

elseif o.plot == 2 %% predictive distribution

    col_=[0.0078 0.6211 0];
    
    t.InnerPosition = [120 120 700 500];
    nexttile(t)
    
    ph(1)=plot(distr.grid.x,squeeze(distr.P_x_G_nu_q(:,:,:,idx_nu,idx_nu)),'k','LineWidth',LW_UB)
    hold on

    % AP
    ph(2)=plot(distr.grid.x,squeeze(distr.P_x_G_nu_q(:,:,:,idx_nu,idx_q)),'Color',color_AP,'LineWidth',LW);
    % ND
    ph(3)=plot(distr.grid.x,squeeze(distr.P_x_G_nu_q(:,:,:,idx_q,idx_nu)),':','Color',color_ND,'LineWidth',LW)
    
    
    xlim([-1 1]*5)
    xticks([-4:2:4])
    yticks([])
    
    box off
    set(gca,'FontSize',FS,'FontName','Calibri Light','layer','top','LineWidth',2)
    exportgraphics(t,'figures\static_predictive_distr.emf','ContentType','vector')

elseif o.plot == 3 %% biases

    t.InnerPosition = [120 120 700 500];
    nexttile(t)
    hold on

    grid_yz = [-fliplr(distr.grid.y) distr.grid.y];
    distr.lP_x_G_yz  = ([fliplr(distr.lP_x_G_y_z(:,:,1)) distr.lP_x_G_y_z(:,:,2)]);
    lP_z_G_x_nu_q    = log(distr.P_x_z_G_nu_q)-log(distr.P_x_G_nu_q);


    beta_=par.beta;
    P_r1_G_x_nu_q   = exp(lP_z_G_x_nu_q(:,:,2,:,:)); P_r1_G_x_nu_q = P_r1_G_x_nu_q.^beta_./(P_r1_G_x_nu_q.^beta_+(1-P_r1_G_x_nu_q).^beta_);
    P_r1_G_x_nu_q(exp(lP_z_G_x_nu_q(:,:,2,:,:))==.5)   = .5;
    P_r1_G_yz_nu_qu = squeeze(sum(P_r1_G_x_nu_q .* exp(distr.lP_x_G_yz),1));


    set(gca,'FontSize',FS,'FontName','Calibri Light','layer','top','LineWidth',2)


    fun = @(x,y_) 1./(1+exp(-(x(1)+x(2)*y_)));
    x0 = [100,-1];


    rp=mean(P_r1_G_yz_nu_qu(80:81,idx_nu,idx_q));
    xAP = log(rp/(1-rp));
    rp=mean(P_r1_G_yz_nu_qu(80:81,idx_q,idx_nu));
    xND = log(rp/(1-rp));
    rp=mean(P_r1_G_yz_nu_qu(80:81,idx_nu,idx_nu));
    xUB = log(rp/(1-rp));

    W=12;
    w=.25;
    h=.025;
    xbase = [[.5 .5]-w [.5 .5]+w];
    fill(xbase+w/2,[0 1 1 0]*xAP(1)+h,color_AP,'EdgeColor','k','LineWidth',W)
    fill(xbase+1-w/2,[0 1 1 0]*xND(1)-h,color_ND,'EdgeColor','k','LineWidth',W)
    fill(xbase+2+w/2,[0 1 1 0]*-xAP(1)-h,'k','EdgeColor',color_AP,'LineWidth',W)
    fill(xbase+3-w/2,[0 1 1 0]*-xND(1)+h,'k','EdgeColor',color_ND,'LineWidth',W)
    plot([0 4],[0 0],'k','LineWidth',3)

    t.InnerPosition = [120 120 300 500];
    axis([0 4 -1 1])
    xticks([])
    exportgraphics(t,'figures\static_psychophysics_biases.emf','ContentType','vector')




elseif o.plot == 4 %% Bias

    t.InnerPosition = [120 120 500 520];
    nexttile(t)
    hold on

    grid_yz = [-fliplr(distr.grid.y) distr.grid.y];
    distr.lP_x_G_yz  = ([fliplr(distr.lP_x_G_y_z(:,:,1)) distr.lP_x_G_y_z(:,:,2)]);
    lP_z_G_x_nu_q    = log(distr.P_x_z_G_nu_q)-log(distr.P_x_G_nu_q);


    %     P_r1_G_x_nu_q   = round(exp(lP_z_G_x_nu_q(:,:,2,:,:)));
    beta_=par.beta;
    P_z_G_x_nu_q  = (exp(lP_z_G_x_nu_q(:,:,2,:,:)));
    P_r1_G_x_nu_q  = P_z_G_x_nu_q.^beta_./(P_z_G_x_nu_q.^beta_+(1-P_z_G_x_nu_q).^beta_);
    P_r1_G_x_nu_q(P_z_G_x_nu_q==.5)   = .5;
    P_r1_G_yz_nu_qu = squeeze(sum(P_r1_G_x_nu_q .* exp(distr.lP_x_G_yz),1));


    rp    = squeeze(mean(P_r1_G_yz_nu_qu(80:81,:,:)));
    bias_ = log(rp./(1-rp));


    J = customcolormap_preset('red-white-blue');
    colormap(flipud(J))
    B = bias_(crop_idx(1):crop_idx(2),crop_idx(1):crop_idx(2));
    pch = pcolor(grid_q,grid_nu,B);
    pch.EdgeColor = 'none';

    view(0,90)
    colorbar('box','off')
    shading interp
    set(gca, 'clim', [-1 1]);


    LHcolor = [1 1 1]*.5;
    Z = avglLH_diff(crop_idx(1):crop_idx(2),crop_idx(1):crop_idx(2));
    contour(grid_q,grid_nu,Z,min(Z(:))+(max(Z(:))-min(Z(:)))*[0:.2:.9 .995],'LineColor',LHcolor(1,:)*.8+[1 1 1]*0,'LineWidth',3)%,'ShowText','on')
    view(0,90)
    scatter3([distr.grid.q(idx_q)],[distr.grid.nu(idx_nu)],[.1],MS,color_AP,'filled','LineWidth',2,'MarkerEdgeColor',[1 1 1])%color_AP*.6)
    scatter3([distr.grid.nu(idx_nu)],[distr.grid.q(idx_q)],[.1],MS,color_ND,'filled','LineWidth',2,'MarkerEdgeColor',[1 1 1])%color_ND*.6)
    scatter3([.5],[.5],[.1],MS,'k','LineWidth',10,'Marker','+')


    axis([crop crop])
    axis square
    xticks([.35 .5 .65])
    yticks([.35 .5 .65])

    set(gca,'FontSize',FS,'FontName','Calibri Light','layer','top','LineWidth',2)
    exportgraphics(t,'figures\static_biases.emf','ContentType','vector')


elseif 5

    crop         = [.4 .7];
    disp(['crop: ',mat2str(crop)])
    [~,crop_idx] = min(abs(distr.grid.q(:)-crop));
    crop_idx     = crop_idx + [-1,1];
    crop_idx = min( crop_idx, length(distr.grid.q)); crop_idx = max( crop_idx, 1);
    grid_q  = squeeze(distr.grid.q(crop_idx(1):crop_idx(2)));
    grid_nu = squeeze(distr.grid.nu(crop_idx(1):crop_idx(2)));

    t.InnerPosition = [120 120 500 520];
    nexttile(t)
    hold on

    grid_yz = [-fliplr(distr.grid.y) distr.grid.y];
    distr.lP_x_G_yz  = ([fliplr(distr.lP_x_G_y_z(:,:,1)) distr.lP_x_G_y_z(:,:,2)]);
    lP_z_G_x_nu_q    = log(distr.P_x_z_G_nu_q)-log(distr.P_x_G_nu_q);

    beta_=par.beta;
    P_r1_G_x_nu_q   = (exp(lP_z_G_x_nu_q(:,:,2,:,:))); P_r1_G_x_nu_q = P_r1_G_x_nu_q.^beta_./(P_r1_G_x_nu_q.^beta_+(1-P_r1_G_x_nu_q).^beta_);
    P_r1_G_x_nu_q(exp(lP_z_G_x_nu_q(:,:,2,:,:))==.5)   = .5;
    P_r1_G_yz_nu_qu = squeeze(sum(P_r1_G_x_nu_q .* exp(distr.lP_x_G_yz),1));

    rp    = squeeze(mean(P_r1_G_yz_nu_qu(80:81,:,:)));
    bias_ = log(rp./(1-rp));

    % J = lines(5)
    J = customcolormap_preset('red-white-blue');
    colormap(flipud(J))
    B = bias_(crop_idx(1):crop_idx(2),crop_idx(1):crop_idx(2));
    pch = pcolor(grid_q,grid_nu,B)
    pch.EdgeColor = 'none';

    view(0,90)
    colorbar('box','off')
    shading interp
%     set(gca, 'clim', [-1 1]);

    scatter3([.5],[.5],[.1],MS,'k','LineWidth',10,'Marker','+')
    plot([1 0],[0 1]+.15,'k','LineWidth',5)
%     plot([1 0],[0 1],'k','LineWidth',5)

    axis([crop crop])
    axis square
    xticks([.35 .5 .65])
    yticks([.35 .5 .65])
%     xticks([])
    set(gca,'FontSize',FS,'FontName','Calibri Light','layer','top','LineWidth',2)
    exportgraphics(t,'figures\demo_bias_1.emf','ContentType','vector')

end
%% LLH comp
return
%%






























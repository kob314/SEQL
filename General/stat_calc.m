function beh_stat = stat_calc(o,data_in)   % behavioural statistics

for ie = o.exp_type_v

    [beh_stat(ie).y_M  beh_stat(ie).dprime_M beh_stat(ie).bias_M] = deal(nan(max(data_in(ie).incl_subj),o.n_bins_y));
    [beh_stat(ie).yz_M beh_stat(ie).r_M beh_stat(ie).r_sim_M] = deal(nan(max(data_in(ie).incl_subj),o.n_bins_y*2));
    [beh_stat(ie).STSE beh_stat(ie).STSE_sim]  = deal(nan(max(data_in(ie).incl_subj),o.n_back));


    for is = data_in(ie).incl_subj
        
        data = data_in(ie).subj(is);

        y_v   = data.y(o.interval);  % stimulus strength
        z_v   = data.z(o.interval);  % trial
        yz_v  = data.yz(o.interval); % signed stimulus strength
        r_v   = data.r(o.interval);  % response

        bins_y    = prctile(y_v,linspace(0,100,o.n_bins_y+1));
        bins_y([1,end]) = [0,1];
        bins_yz  = [-fliplr(bins_y) bins_y(2:end)];

        cats_y = discretize(y_v,bins_y);
        cats_yz = discretize(yz_v,bins_yz);

        %% psychometric curve
        yz_ = accumarray(cats_yz(:),yz_v(:),[o.n_bins_y*2 1],@nanmean);
        r_  = accumarray(cats_yz(:),r_v(:),[o.n_bins_y*2 1],@nanmean);

        beh_stat(ie).yz_M(is,:)  = yz_;
        beh_stat(ie).r_M(is,:)   = r_;


        %% SDT dprime & bias
        y_  = accumarray(cats_y',y_v',[o.n_bins_y 1],@nanmean);
        beh_stat(ie).y_M(is,:) = accumarray(cats_y(:),y_v(:),[o.n_bins_y 1],@nanmean);

        for ic = 1:o.n_bins_y
            [beh_stat(ie).dprime_M(is,ic) beh_stat(ie).bias_M(is,ic)] = dp_bias_calc(z_v(cats_y==ic), r_v(cats_y==ic));
        end


        %% short-term serial effect (STSE)

        past_r_m = data.past_r(:,o.interval);
        past_z_m = data.past_z(:,o.interval);


        A = nansum(r_v.*past_r_m    .*past_z_m,2)    ./nansum(~isnan(r_v).*past_r_m.*past_z_m,2);
        B = nansum(r_v.*past_r_m    .*(1-past_z_m),2)./nansum(~isnan(r_v).*past_r_m.*(1-past_z_m),2);
        C = nansum(r_v.*(1-past_r_m).*past_z_m,2)    ./nansum(~isnan(r_v).*(1-past_r_m).*past_z_m,2);
        D = nansum(r_v.*(1-past_r_m).*(1-past_z_m),2)./nansum(~isnan(r_v).*(1-past_r_m).*(1-past_z_m),2);

        beh_stat(ie).STSE(is,:) = (nanmean([A,B],2)-nanmean([C,D],2))/2+1/2;

        A = nansum(z_v.*past_z_m    .*past_r_m,2)    ./nansum(~isnan(z_v).*past_z_m.*past_r_m,2);
        B = nansum(z_v.*past_z_m    .*(1-past_r_m),2)./nansum(~isnan(z_v).*past_z_m.*(1-past_r_m),2);
        C = nansum(z_v.*(1-past_z_m).*past_r_m,2)    ./nansum(~isnan(z_v).*(1-past_z_m).*past_r_m,2);
        D = nansum(z_v.*(1-past_z_m).*(1-past_r_m),2)./nansum(~isnan(z_v).*(1-past_z_m).*(1-past_r_m),2);

        beh_stat(ie).STSE_z(is,:) = (nanmean([A,B],2)-nanmean([C,D],2))/2+1/2;


    end
end

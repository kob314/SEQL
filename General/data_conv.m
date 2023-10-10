function data_out = data_conv(data_in,o)


for ie = o.exp_type_v

    session_length = length(data_in.test__.exp_num(ie).subj_num);
    inc_v = nan(1,session_length);
    inc0_v = nan(1,session_length);
    inc1_v = nan(1,session_length);
    for is = 1:session_length

        data_out(ie).subj(is).z       = data_in.test__.exp_num(ie).subj_num(is).trial_v;
        data_out(ie).subj(is).r       = data_in.test__.exp_num(ie).subj_num(is).response_v;
%         data_out(ie).subj(is).r       = data_out(ie).subj(is).z;
%         data_out(ie).subj(is).r       = data_out(ie).subj(is).r( randperm( length(data_out(ie).subj(is).r) ) )
%         aux = rand(size(data_out(ie).subj(is).z)) < .8;
%         data_out(ie).subj(is).r       = data_out(ie).subj(is).z .* aux +(1-aux) .* rand(size(data_out(ie).subj(is).z))
        data_out(ie).subj(is).correct = 1-abs(data_out(ie).subj(is).z - data_out(ie).subj(is).r);
        data_out(ie).subj(is).rt      = data_in.test__.exp_num(ie).subj_num(is).reac_t_v;

        if isfield(data_in.train2.exp_num(ie).subj_num(is),'trial_v')
            data_out(ie).last_jump_type(is) = mean( data_in.train2.exp_num(ie).subj_num(is).trial_v(70:90)) > .5; 
        else
            data_out(ie).last_jump_type(is) = nan;
        end

        for i=1:2

            switch i
                case 1
                    aux_v = 2*(data_out(ie).subj(is).r-.5);
                case 2
                    aux_v = 2*(data_out(ie).subj(is).z-.5);
            end
            aux=repmat([aux_v(:); 0],o.n_back+1,1);
            aux=tril(reshape( aux(1:(length(aux_v)*(o.n_back+1))) ,length(aux_v),[]));
            aux = aux(:,2:end);
            aux(aux==0) = NaN;
            aux = (aux+1)/2;
            switch i
                case 1
                    data_out(ie).subj(is).past_r = aux';
                case 2
                    data_out(ie).subj(is).past_z = aux';
            end

        end




        %% noise rescale
        gamma_v = data_in.test__.exp_num(ie).subj_num(is).gamma_v;
        if max(gamma_v) > 100
            gamma_v = gamma_v/10;
        end
        gamma_v(gamma_v<1) = 1;  %% fitting is such
%         gamma_v = gamma_v-1;
%         data_out(ie).subj(is).y = 1./(1+gamma_v);
        data_out(ie).subj(is).y = 1./gamma_v;

        data_out(ie).subj(is).y_min_max = [min(data_out(ie).subj(is).y) max(data_out(ie).subj(is).y)];

        data_out(ie).subj(is).yz = 2*(data_out(ie).subj(is).z-.5).*data_out(ie).subj(is).y;

        y_bins = [0    0.0555    0.0682    0.0888    0.1277    0.2267    1];
        y_bins(1) = 0;
        psych(ie).ys_bins(is,:)= y_bins;

        yz_bins = [-fliplr(y_bins) y_bins(2:end)];
        data_out(ie).subj(is).yz_cats = discretize(data_out(ie).subj(is).yz,yz_bins);

        %% filtering subjects
        X = data_out(ie).subj(is).yz(~isnan(data_out(ie).subj(is).r));
        Y = data_out(ie).subj(is).r(~isnan(data_out(ie).subj(is).r));
        [rho,plev]=corr(X(:),Y(:));
        data_out(ie).subj(is).correlation_v = [rho,plev];
        
        X0 = data_out(ie).subj(is).y(  ~isnan(data_out(ie).subj(is).r)  &  data_out(ie).subj(is).z == 0  );
        X1 = data_out(ie).subj(is).y(  ~isnan(data_out(ie).subj(is).r)  &  data_out(ie).subj(is).z == 1  );
        Y0 = data_out(ie).subj(is).correct(  ~isnan(data_out(ie).subj(is).r)  &  data_out(ie).subj(is).z == 0  );
        Y1 = data_out(ie).subj(is).correct(  ~isnan(data_out(ie).subj(is).r)  &  data_out(ie).subj(is).z == 1  );

        if ~isempty(X0) & ~isempty(Y0) & ~isempty(X1) & ~isempty(Y1)

            [rho,plev]=corr(X0(:),Y0(:));
            data_out(ie).subj(is).correlation0_v = [rho,plev];

            [rho,plev]=corr(X1(:),Y1(:));
            data_out(ie).subj(is).correlation1_v = [rho,plev];


            inc_v(is) = (data_out(ie).subj(is).correlation_v(1)) > 0 & (data_out(ie).subj(is).correlation_v(2) < 0.05);
            inc0_v(is) = data_out(ie).subj(is).correlation0_v(1) > 0;
            inc1_v(is) = data_out(ie).subj(is).correlation1_v(1) > 0;
%             inc0_v(is) = (data_out(ie).subj(is).correlation0_v(1) > 0) & (data_out(ie).subj(is).correlation0_v(2) < 0.1);
%             inc1_v(is) = (data_out(ie).subj(is).correlation1_v(1) > 0) & (data_out(ie).subj(is).correlation1_v(2) < 0.1);

           
            data_out(ie).subj(is).obj_pair = data_in.exp_num(ie).subj_num(is).obj_pair;
        else
            inc_v  = data_out(ie).subj(is).z*0;
            inc0_v = inc_v;
            inc1_v = inc_v;
            data_out(ie).subj(is).obj_pair = [];
        end
     


    end

    data_out(ie).incl_subj = find(inc_v & inc0_v & inc1_v);
    data_out(ie).excl_subj = find(~(inc_v & inc0_v & inc1_v));
%     data_in
    data_out(ie).train.AP       = data_in.train1.exp_num(ie).p;
    data_out(ie).test.AP        = data_in.test__.exp_num(ie).p;



end


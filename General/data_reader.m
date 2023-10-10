function data_ = data_reader

% load('Simdata_.mat')
% TestData_ = TestData;
load('AllDat_ext.mat')
% numsubj = [21 0 18 19];
% for ii = [1,3,4]
%     try
%         TestData(end+1,1:size(TestData_(ii,:),2)) = TestData_(ii,:);
%         NumSubj =[NumSubj numsubj(ii)];
%     end
% end


%% load experiment data
% variable: TrainData, TestData
% 1. row: Answer 1 or 2
% 2. row: Presented object 1 or 2
% 3. row: Reaction time
% 4. row: Noise on trial (noise*10)


% original {"50-65", "65-65", "75-65", "50-65Grad", "50-65Catch", "50-50Pilot", "50-6550Decr", "MidShift75", "75-50", "Guess5065", "50-7550Grad",'Asym 50-50','Sym 50-65','Sym 50-75-50','Vol 50-50-65', 'Vol 75-75-50', 'Stat 75-75-50'}
% new list {'50-65','75-50','75-65','65-65','sim','50-65 ramp','50-65-50 ramp','50-75-50 ramp','50-75 mid shift','50-65Catch','50-65Guess','Asym 50-50','Sym 50-65','Sym 50-75-50','Vol 50-50-65','Vol 75-75-50', 'Stat 75-75-50'};

% ExpNames{end+1} = 'sim 50-65';
% ExpNames{end+1} = 'sim 50-65';
% ExpNames{end+1} = 'sim 75-65';
% ExpNames{end+1} = 'sim 65-65';

% ---------------- 1 2 3 4      5           6 7 8  9 10 11 12 13 14 15 16
% 17
% ExpNames
% ie_list    = [1 9 3 2 14 4 7 11 8 5   10 6 12 13 14 15 16];
ie_list      = [1 9 3 2 16 4 7 11 8 5   10 6 12 13 14 15 16 17 18 1 2];
exp_names       = ExpNames(ie_list);
for ie = 1:length(ie_list)
    exp_n = ie_list(ie);
    
    switch ie
        case 1
            % exp_name = 'exp4Unbal.mat';
            % [30:35,95:102,104:108,130,131,132,150:154]
            data_.exp_num(ie).N_subj = NumSubj(exp_n);
            data_.train1.exp_num(ie).p = .5;
            data_.test__.exp_num(ie).p    = .65;
        case 2
            % exp_name = 'exp4TrainUnbaltest__B.mat';
            % subj_v = setdiff([140:155,157:159,161],excl);  % without [157:159,161]
            data_.exp_num(ie).N_subj = NumSubj(exp_n);
            data_.train1.exp_num(ie).p = .75;
            data_.test__.exp_num(ie).p    = .5;
        case 3
            % exp_name = 'exp4UnbalTrained.mat';
            % subj_v = setdiff([112:129],excl);  % previously 111:123
            data_.exp_num(ie).N_subj = NumSubj(exp_n);
            data_.train1.exp_num(ie).p = .75;
            data_.test__.exp_num(ie).p    = .65;
        case 4
            % exp_name = 'exp46565.mat';
            % subj_v = setdiff([250:256,258:269],excl);
            data_.exp_num(ie).N_subj = NumSubj(exp_n);
            data_.train1.exp_num(ie).p = .65;
            data_.test__.exp_num(ie).p    = .65;
        case 5
            data_.exp_num(ie).N_subj = NumSubj(exp_n);
            data_.train1.exp_num(ie).p = .5;
            data_.test__.exp_num(ie).p    = .65;
        case 6  % ramping 50->65
            %exp_name = 'ExpIVGrad.mat';
            %subj_v = setdiff([301:319],excl);
            data_.exp_num(ie).N_subj = NumSubj(exp_n);
            data_.train1.exp_num(ie).p = .5;
            data_.test__.exp_num(ie).p    = .65;
        case 7 % ramping 50-65->50
            % exp_name = 'expIVGradDecr.mat';
            % subj_v = setdiff([401:414,416,417,419:422],excl); %setdiff([401:403,405:408,411,419],excl);
            data_.exp_num(ie).N_subj = NumSubj(exp_n);
            data_.train1.exp_num(ie).p = .5;
            data_.test__.exp_num(ie).p    = .5;
        case 8 % ramping 50-75-50>
            data_.exp_num(ie).N_subj = NumSubj(exp_n);
            data_.train1.exp_num(ie).p = .5;
            data_.test__.exp_num(ie).p    = .5;
        case 9 % midshift
            data_.exp_num(ie).N_subj = NumSubj(exp_n);
            data_.train1.exp_num(ie).p = .5;
            data_.test__.exp_num(ie).p    = .75;
        case 10
            % exp_name = 'exp4Catch.mat';
            % subj_v = setdiff([200:202,204:218],excl);
            data_.exp_num(ie).N_subj = NumSubj(exp_n);
            data_.train1.exp_num(ie).p = .5;
            data_.test__.exp_num(ie).p    = .65;
        case 11 % midshift
            % exp_name = 'expIV5065MixGuess.mat';
            % subj_v = setdiff([101:118],excl);
            data_.exp_num(ie).N_subj = NumSubj(exp_n);
            data_.train1.exp_num(ie).p = .5;
            data_.test__.exp_num(ie).p    = .65;
        case 12 % fifty-fifty
            data_.exp_num(ie).N_subj = NumSubj(exp_n);
            data_.train1.exp_num(ie).p = .5;
            data_.test__.exp_num(ie).p    = .5;
        case 13 % Asy fifty-fifty
            data_.exp_num(ie).N_subj = NumSubj(exp_n);
            data_.train1.exp_num(ie).p = .5;
            data_.test__.exp_num(ie).p    = .5;
        case 14 % Sym 50-65
            data_.exp_num(ie).N_subj = NumSubj(exp_n);
            data_.train1.exp_num(ie).p = .5;
            data_.test__.exp_num(ie).p    = .65;
        case 15 % Sym 50-75-50
            data_.exp_num(ie).N_subj = NumSubj(exp_n);
            data_.train1.exp_num(ie).p = .75;
            data_.test__.exp_num(ie).p    = .5;
        case 16 % Vol 50-50-65
            data_.exp_num(ie).N_subj = NumSubj(exp_n);
            data_.train1.exp_num(ie).p = .5;
            data_.test__.exp_num(ie).p    = .65;
        case 17 % Vol 50-50-65b
            data_.exp_num(ie).N_subj = NumSubj(exp_n);
            data_.train1.exp_num(ie).p = .5;
            data_.test__.exp_num(ie).p    = .65;
        case 18 % Vol 75-75-50
            data_.exp_num(ie).N_subj = NumSubj(exp_n);
            data_.train1.exp_num(ie).p = .75;
            data_.test__.exp_num(ie).p    = .5;
        case 19 % Stat 75-75-50
            data_.exp_num(ie).N_subj = NumSubj(exp_n);
            data_.train1.exp_num(ie).p = .75;
            data_.test__.exp_num(ie).p    = .5;
        case 19 % Stat 75-75-50
            data_.exp_num(ie).N_subj = NumSubj(exp_n);
            data_.train1.exp_num(ie).p = .75;
            data_.test__.exp_num(ie).p    = .5;
        case 20
            load('AllDat_ext_tilt.mat')
            data_.exp_num(ie).N_subj = NumSubj(1);
            data_.train1.exp_num(ie).p = .5;
            data_.test__.exp_num(ie).p    = .65;
            exp_names{20} = '50-65_tilt';
        case 21
            load('AllDat_ext_tilt.mat')
            data_.exp_num(ie).N_subj = NumSubj(2);
            data_.train1.exp_num(ie).p = .65;
            data_.test__.exp_num(ie).p    = .65;
            exp_names{21} = '65-65_tilt';
%         case 15 % Simm 50-65
%             data_.exp_num(ie).N_subj = NumSubj(exp_n);
%             data_.train1.exp_num(ie).p = .5;
%             data_.test__.exp_num(ie).p    = .65;
%         case 16 % Simm 75-65
%             data_.exp_num(ie).N_subj = NumSubj(exp_n);
%             data_.train1.exp_num(ie).p = .75;
%             data_.test__.exp_num(ie).p    = .65;
%         case 17 % Simm 65-65
%             data_.exp_num(ie).N_subj = NumSubj(exp_n);
%             data_.train1.exp_num(ie).p = .65;
%             data_.test__.exp_num(ie).p    = .65;
    end
    
%     NumSubj
    for is = 1:NumSubj(exp_n)
        
        
        
        if ie == 2
            if mean( TrainData{exp_n,is}(2,:)-1 ) > .5
                data_.train1.exp_num(ie).subj_num(is).response_v = TrainData{exp_n,is}(1,:)-1;
                data_.test__.exp_num(ie).subj_num(is).response_v  = TestData{exp_n,is}(1,:)-1;
                data_.train1.exp_num(ie).subj_num(is).trial_v    = TrainData{exp_n,is}(2,:)-1;
                data_.test__.exp_num(ie).subj_num(is).trial_v     = TestData{exp_n,is}(2,:)-1;
                data_.exp_num(ie).subj_num(is).obj_pair = Objects{exp_n,is};
            else
                data_.train1.exp_num(ie).subj_num(is).response_v  = 2-TrainData{exp_n,is}(1,:);
                data_.test__.exp_num(ie).subj_num(is).response_v   = 2-TestData{exp_n,is}(1,:);
                data_.train1.exp_num(ie).subj_num(is).trial_v     = 2-TrainData{exp_n,is}(2,:);
                data_.test__.exp_num(ie).subj_num(is).trial_v      = 2-TestData{exp_n,is}(2,:);
                data_.exp_num(ie).subj_num(is).obj_pair = fliplr(Objects{exp_n,is});
            end
        elseif ismember(ie,[7,8])
            if mean( TestData{exp_n,is}(2,1:60)-1 ) > .5
                data_.train1.exp_num(ie).subj_num(is).response_v = TrainData{exp_n,is}(1,:)-1;
                data_.test__.exp_num(ie).subj_num(is).response_v  = TestData{exp_n,is}(1,:)-1;
                data_.train1.exp_num(ie).subj_num(is).trial_v    = TrainData{exp_n,is}(2,:)-1;
                data_.test__.exp_num(ie).subj_num(is).trial_v     = TestData{exp_n,is}(2,:)-1;
                data_.exp_num(ie).subj_num(is).obj_pair = Objects{exp_n,is};
            else
                data_.train1.exp_num(ie).subj_num(is).response_v  = 2-TrainData{exp_n,is}(1,:);
                data_.test__.exp_num(ie).subj_num(is).response_v   = 2-TestData{exp_n,is}(1,:);
                data_.train1.exp_num(ie).subj_num(is).trial_v     = 2-TrainData{exp_n,is}(2,:);
                data_.test__.exp_num(ie).subj_num(is).trial_v      = 2-TestData{exp_n,is}(2,:);
                data_.exp_num(ie).subj_num(is).obj_pair = fliplr(Objects{exp_n,is});
            end
        elseif ismember(ie,[5])%,15,16,17])        
%             try
%                 'try'
                data_.test__.exp_num(ie).subj_num(is).response_v  = TestData{exp_n,is}(1,:);
                data_.test__.exp_num(ie).subj_num(is).trial_v     = TestData{exp_n,is}(2,:);
                data_.test__.exp_num(ie).subj_num(is).reac_t_v    = TestData{exp_n,is}(3,:);
                data_.test__.exp_num(ie).subj_num(is).gamma_v     = TestData{exp_n,is}(4,:);
                data_.test__.exp_num(ie).subj_num(is).min_max = [min(data_.test__.exp_num(ie).subj_num(is).gamma_v),max(data_.test__.exp_num(ie).subj_num(is).gamma_v)];
                data_.test__.exp_num(ie).subj_num(is).catch_idx   = [];
%             end
        elseif ismember(ie,[13:19])
%             if ie == 13
%                 crit_ = mean( TrainNData{exp_n,is}(4,TrainNData{exp_n,is}(2,:)==2 ) ) < mean( TrainNData{exp_n,is}(4,TrainNData{exp_n,is}(2,:)==1) );
%             elseif ie == 14
% %                 crit_ = mean( TestData{exp_n,is}(2,:)-1 ) > .5 ;
%                 crit_ = mean( TrainNData{exp_n,is}(2,:)-1 ) > .5 ;
%             end
            switch ie
                case 13
                    crit_ = mean( TrainNData{exp_n,is}(4,TrainNData{exp_n,is}(2,:)==2 ) ) < mean( TrainNData{exp_n,is}(4,TrainNData{exp_n,is}(2,:)==1) );
                case 14
                    crit_ = mean( TestData{exp_n,is}(2,:)-1 ) > .5 ;
                case 15
                    crit_ = mean( TrainNData{exp_n,is}(2,:)-1 ) > .5 ;
                case 16
                    crit_ = mean( TestData{exp_n,is}(2,:)-1 ) > .5 ;
                case 17
                    crit_ = mean( TestData{exp_n,is}(2,:)-1 ) > .5 ;
                case 18
                    crit_ = mean( TrainNData{exp_n,is}(2,:)-1 ) > .5 ;
                case 19
                    crit_ = mean( TrainNData{exp_n,is}(2,:)-1 ) > .5 ;
            end
            
            if crit_
                data_.train1.exp_num(ie).subj_num(is).response_v = TrainData{exp_n,is}(1,:)-1;
                data_.train2.exp_num(ie).subj_num(is).response_v = TrainNData{exp_n,is}(1,:)-1;
                data_.test__.exp_num(ie).subj_num(is).response_v  = TestData{exp_n,is}(1,:)-1;
                data_.train1.exp_num(ie).subj_num(is).trial_v    = TrainData{exp_n,is}(2,:)-1;
                data_.train2.exp_num(ie).subj_num(is).trial_v    = TrainNData{exp_n,is}(2,:)-1;
                data_.test__.exp_num(ie).subj_num(is).trial_v     = TestData{exp_n,is}(2,:)-1;
                data_.exp_num(ie).freq(is) = 2;
                data_.exp_num(ie).subj_num(is).obj_pair = Objects{exp_n,is};
                
            else
                data_.train1.exp_num(ie).subj_num(is).response_v  = 2-TrainData{exp_n,is}(1,:);
                data_.train2.exp_num(ie).subj_num(is).response_v  = 2-TrainNData{exp_n,is}(1,:);
                data_.test__.exp_num(ie).subj_num(is).response_v   = 2-TestData{exp_n,is}(1,:);
                data_.train1.exp_num(ie).subj_num(is).trial_v     = 2-TrainData{exp_n,is}(2,:);
                data_.train2.exp_num(ie).subj_num(is).trial_v     = 2-TrainNData{exp_n,is}(2,:);
                data_.test__.exp_num(ie).subj_num(is).trial_v      = 2-TestData{exp_n,is}(2,:);
                data_.exp_num(ie).freq(is) = 1;
                data_.exp_num(ie).subj_num(is).obj_pair = fliplr(Objects{exp_n,is});
            end
            data_.train1.exp_num(ie).subj_num(is).response_v(data_.train1.exp_num(ie).subj_num(is).response_v==-1)=nan;
            data_.train1.exp_num(ie).subj_num(is).response_v(data_.train1.exp_num(ie).subj_num(is).response_v==2)=nan;
            data_.train2.exp_num(ie).subj_num(is).response_v(data_.train2.exp_num(ie).subj_num(is).response_v==-1)=nan;
            data_.train2.exp_num(ie).subj_num(is).response_v(data_.train2.exp_num(ie).subj_num(is).response_v==2)=nan;
        elseif ismember(ie,[20:21])
            if mean( TestData{exp_n,is}(2,:)-1 ) > .5
                data_.test__.exp_num(ie).subj_num(is).response_v  = TestData{exp_n,is}(1,:)-1;
                data_.test__.exp_num(ie).subj_num(is).trial_v     = TestData{exp_n,is}(2,:)-1;
            else
                data_.test__.exp_num(ie).subj_num(is).response_v   = 2-TestData{exp_n,is}(1,:);
                data_.test__.exp_num(ie).subj_num(is).trial_v      = 2-TestData{exp_n,is}(2,:);
            end
        else
            if mean( TestData{exp_n,is}(2,:)-1 ) > .5
                data_.train1.exp_num(ie).subj_num(is).response_v = TrainData{exp_n,is}(1,:)-1;
                data_.test__.exp_num(ie).subj_num(is).response_v  = TestData{exp_n,is}(1,:)-1;
                data_.train1.exp_num(ie).subj_num(is).trial_v    = TrainData{exp_n,is}(2,:)-1;
                data_.test__.exp_num(ie).subj_num(is).trial_v     = TestData{exp_n,is}(2,:)-1;
                data_.exp_num(ie).subj_num(is).obj_pair = Objects{exp_n,is};
            else
                data_.train1.exp_num(ie).subj_num(is).response_v  = 2-TrainData{exp_n,is}(1,:);
                data_.test__.exp_num(ie).subj_num(is).response_v   = 2-TestData{exp_n,is}(1,:);
                data_.train1.exp_num(ie).subj_num(is).trial_v     = 2-TrainData{exp_n,is}(2,:);
                data_.test__.exp_num(ie).subj_num(is).trial_v      = 2-TestData{exp_n,is}(2,:);
                data_.exp_num(ie).subj_num(is).obj_pair = fliplr(Objects{exp_n,is});
            end
        end
        
        
        if ~ismember(ie,[5])%,15,16,17])
            
            if ismember(ie,[20:21])
                data_.test__.exp_num(ie).subj_num(is).reac_t_v   = TestData{exp_n,is}(3,:);
                data_.test__.exp_num(ie).subj_num(is).gamma_v   = TestData{exp_n,is}(4,:);
            else
                
                data_.train1.exp_num(ie).subj_num(is).reac_t_v  = TrainData{exp_n,is}(3,:);
                data_.test__.exp_num(ie).subj_num(is).reac_t_v   = TestData{exp_n,is}(3,:);
                
                data_.train1.exp_num(ie).subj_num(is).gamma_v  = TrainData{exp_n,is}(4,:);
                data_.test__.exp_num(ie).subj_num(is).gamma_v   = TestData{exp_n,is}(4,:);
                
                data_.train1.exp_num(ie).subj_num(is).min_max = [min(data_.train1.exp_num(ie).subj_num(is).gamma_v),max(data_.train1.exp_num(ie).subj_num(is).gamma_v)];
                data_.test__.exp_num(ie).subj_num(is).min_max  = [min(data_.test__.exp_num(ie).subj_num(is).gamma_v),max(data_.test__.exp_num(ie).subj_num(is).gamma_v)];
                
                
                if ismember(ie,[13:19])
                    data_.train2.exp_num(ie).subj_num(is).reac_t_v  = TrainNData{exp_n,is}(3,:);
                    data_.train2.exp_num(ie).subj_num(is).gamma_v  = TrainNData{exp_n,is}(4,:);
                    data_.train2.exp_num(ie).subj_num(is).min_max = [min(data_.train2.exp_num(ie).subj_num(is).gamma_v),max(data_.train2.exp_num(ie).subj_num(is).gamma_v)];
                    %                 data_.exp_num(ie).objects{is} = Objects{exp_n,is};
                end
                
                if ie == 10
                    data_.train1.exp_num(ie).subj_num(is).catch_idx = [];
                    data_.test__.exp_num(ie).subj_num(is).catch_idx = TestData{exp_n,is}(2,:) == 0;
                elseif ie == 11
                    data_.train1.exp_num(ie).subj_num(is).catch_idx = [];
                    data_.test__.exp_num(ie).subj_num(is).catch_idx = TestData{exp_n,is}(5,:);
                    %                 TestData{exp_n,is}(5,:)
                else
                    data_.train1.exp_num(ie).subj_num(is).catch_idx = [];
                    data_.train2.exp_num(ie).subj_num(is).catch_idx = [];
                    data_.test__.exp_num(ie).subj_num(is).catch_idx = [];
                end
                
            end
            
        end
        
        
    end
end


% 201028 jyh
% make a file for each ROI that contains permuted color selectivity values of all subjects and all cond.
% should fed into a R code later

clear all;


%% Experiment info. & directory

exp='v3';
% ROIs = {'V2_3dg_fmasked_q.001','V3_3dg_fmasked_q.001', 'V4v_fmasked_q.001' }; %, 'V2_3dg_fmasked_q.001','V3_3dg_fmasked_q.001', 'V4v_fmasked_q.001'};
ROIs = {'LGN_hk'}; %, 'V2_3dg_fmasked_q.001','V3_3dg_fmasked_q.001', 'V4v_fmasked_q.001'};
% ROIs = {'V1_3dg_fmasked_q.001'};

nChan = 8; % num. of channels
nColor = 8; % num. of colors used in the experiment
nAttCond = 2; % attention, inattention
nCond = nColor*nAttCond; % num. of condition (8 colors X in vs. out)
nColCond = 2; % cardinal, intercardinal
nLMSCond = 2; % L-M, S

carInx = ([1:(nColor/2)]-1)*2+1; %1, 3, 5, 7
intInx = ([1:(nColor/2)])*2; %2,4,6,8
LMInx = [1,5];
SInx = [3,7];

nBin = 1;

baseDir = sprintf('/Volumes/Duri/data/Color%s', exp);
addpath('/Users/auna/Script/');

fmDir = 'Img_data/forwardmodel'; %
tunDir = fullfile(fmDir, 'sc_dt_hp_am');  %where tuning value txt is located
permDir = fullfile(fmDir, 'perm_sc_dt_hp_am');
statDir = fullfile(baseDir, 'perm_fidelity');
statFileName = 'fidelity_2TRlag_neutralModel.csv';

nPerm = 1000;

%% Subject info.
% %v3
SN = {};
SN{end+1} = '01'; % PSY
SN{end+1} = '03'; % KIS
SN{end+1} = '04'; % JBH
SN{end+1} = '05'; % HJH %%% check alignment for sub4,5 - 6 deg mask 181210
SN{end+1} = '06'; % KB
SN{end+1} = '07'; % SHY
SN{end+1} = '08'; % LSY
SN{end+1} = '09'; % HJH2
SN{end+1} = '10'; % CES
SN{end+1} = '11'; % YYH
SN{end+1} = '12'; % LHB
SN{end+1} = '13'; % JHY
SN{end+1} = '14'; % KYJ


%% Parameters for plots
lineWidth = 3;
lineColors = {[1 0.6 0.784],...
    [1 0 0.6],...
    [0.855 0.702 1],...
    [0.18 0.40 0.73],...
    [0 1 1],...
    [0 0.498 0],...
    [0.4 0.8 0],...
    [0.878 0.537 0.098]};

Hnames= {'All colors', 'Cardinal colors only'};
CInames = {'Cardinal', 'Int.cardinal'};
LMSnames = {'L-M', 'S'};
Attnames = {'In', 'Out'};
car_LMnames = {'L-M color 1', 'L-M color 5'};
car_Snames = {'S color 3', 'S color 7'};

CIColors = {lineColors{2}, lineColors{4}; lineColors{2}, lineColors{4}};
LMSColors = {lineColors{2}, lineColors{4}; lineColors{2}, lineColors{4}};
car_LMColors = {lineColors{LMInx(1)},...
    lineColors{LMInx(1)}*0.6;lineColors{LMInx(2)}, lineColors{LMInx(2)}*0.6};
car_SColors = {lineColors{SInx(1)},...
    lineColors{SInx(1)}*0.6;lineColors{SInx(2)}, lineColors{SInx(2)}*0.6};
opacity = 0.3;
%% Selectivity index
Methods={'Information fidelity'};
sigtxt = {'not significant', 'significant'};

%% permutation array
perm_colSel = cell(nPerm, length(Methods), length(Hnames), 2,2);
perm_carLMS_colSel = cell(nPerm, length(Methods), length(Hnames), 2,2);
perm_ttfit = zeros(nPerm, length(Methods), length(Hnames), 2, 2, 3); %last dim: p,t, & sig value
perm_carLMS_ttfit = zeros(nPerm, length(Methods), length(Hnames), 2, 2, 3); %last dim: p,t, & sig value
perm_tval = cell(length(ROIs), length(Methods), length(Hnames), 2,2);
perm_carLMS_tval = cell(length(ROIs), length(Methods), length(Hnames), 2,2);

perm_rmfit = zeros(nPerm, length(ROIs), length(Methods), length(Hnames), 3);  %main1,2,int1, %last dim: p,f, & sig value
perm_carLMS_rmfit = zeros(nPerm, length(ROIs), length(Methods), length(Hnames),3);
perm_fval = cell(length(ROIs), length(Methods), length(Hnames)); %main1, main2, int
perm_carLMS_fval = cell(length(ROIs), length(Methods), length(Hnames), 3);
perm_ANOVA = cell(length(Methods), length(Hnames));

multcomp_test = zeros(nPerm, length(ROIs), length(Methods), length(Hnames), 4);
carLMS_multcomp_test = zeros(nPerm, length(ROIs), length(Methods), length(Hnames), 4);

factorNames = {'C vs. I ', 'Attention'; 'L-M vs. S', 'Attention'};
carLMS_factorNames = {'L-M1 vs. L-M5', 'Attention'; 'S3 vs. S7', 'Attention'};


%% real_data array
real_tval = cell(length(ROIs), length(Methods), length(Hnames), 2,2);
real_carLMS_tval = cell(length(ROIs), length(Methods), length(Hnames), 2,2);


rmfit = zeros(length(ROIs), length(Methods), length(Hnames),3);
carLMS_rmfit = zeros(length(ROIs), length(Methods), length(Hnames),3);

ANOVA = cell(length(Methods), length(Hnames));

multcomp_test = zeros(length(ROIs), length(Methods), length(Hnames), 4);
carLMS_multcomp_test = zeros(length(ROIs), length(Methods), length(Hnames), 4);


%% load perm data
for xROI = 1:length(ROIs)
        fprintf('.....ROI: %s permutation start.....\n\n', ROIs{xROI})
    for xPerm = 1:nPerm
                if mod(xPerm,100) == 0; fprintf('ROI: %s testing no.%d\n', ROIs{xROI}, xPerm); end
        for xSN = 1:length(SN)
            fileName = sprintf('perm%d_%s_tuning_%s2TR_used4TRs_neutralModel_shift.txt', xPerm, SN{xSN}, ROIs{xROI});
            TT{xSN} = load(fullfile(baseDir, SN{xSN}, permDir, [ROIs{xROI} '_2TRlag_neutralModel'], fileName));
            
            % zscore across channel
            for rr = 1:nCond % cc: row, color x attention
                zTT{xSN}(rr,1:nChan) = zscore(TT{xSN}(rr, :));
            end
        end
        
        
        tTT = [];
        for xSN = 1:length(SN)
            for rr = 1:nCond
                tTT(xSN, 1:nChan, rr) = zTT{xSN}(rr,:);
            end
        end
        
        
        cTT = [];
        cTT(:,1:nChan,1,1) = mean(tTT(:,:,carInx), 3); %cardinal, in
        cTT(:,1:nChan,2,1) = mean(tTT(:,:,intInx), 3); %intercardinal, in
        cTT(:,1:nChan,1,2) = mean(tTT(:,:,carInx+nColor), 3); %cardinal, out
        cTT(:,1:nChan,2,2) = mean(tTT(:,:,intInx+nColor), 3); %intercardinal, out
        car_cTT = [];
        car_cTT(:,1:nChan,1,1) = mean(tTT(:,:,LMInx), 3); %L-M color, in
        car_cTT(:,1:nChan,2,1) = mean(tTT(:,:,SInx), 3); %S color, in
        car_cTT(:,1:nChan,1,2) = mean(tTT(:,:,LMInx+nColor), 3); %L-M color, out
        car_cTT(:,1:nChan,2,2) = mean(tTT(:,:,SInx+nColor), 3); %S color, out
        carLM_cTT = [];
        carLM_cTT(:,1:nChan,1,1) = mean(tTT(:,:,LMInx(1)), 3); %L-M color1, in
        carLM_cTT(:,1:nChan,2,1) = mean(tTT(:,:,LMInx(2)), 3); %L-M color2, in
        carLM_cTT(:,1:nChan,1,2) = mean(tTT(:,:,LMInx(1)+nColor), 3); %L-M color1, out
        carLM_cTT(:,1:nChan,2,2) = mean(tTT(:,:,LMInx(2)+nColor), 3); %L-M color2, out
        carS_cTT = [];
        carS_cTT(:,1:nChan,1,1) = mean(tTT(:,:,SInx(1)), 3); %L-M color1, in
        carS_cTT(:,1:nChan,2,1) = mean(tTT(:,:,SInx(2)), 3); %L-M color2, in
        carS_cTT(:,1:nChan,1,2) = mean(tTT(:,:,SInx(1)+nColor), 3); %L-M color1, out
        carS_cTT(:,1:nChan,2,2) = mean(tTT(:,:,SInx(2)+nColor), 3); %L-M color2, out
        
        
        % Add one more channel to the end to make it odd number
        % Add the last channel value as channel 9
        pTT = []; car_pTT = []; carLM_pTT = []; carS_pTT = [];
        
        pTT(:,1+1:nChan+1,:,:) = cTT(:,1:nChan,:,:);
        pTT(:,1,:,:) = pTT(:,nChan+1,:,:);
        
        car_pTT(:,1+1:nChan+1,:,:) = car_cTT(:,1:nChan,:,:);
        car_pTT(:,1,:,:) = car_pTT(:,nChan+1,:,:);
        
        carLM_pTT(:,1+1:nChan+1,:,:) = carLM_cTT(:,1:nChan,:,:);
        carLM_pTT(:,1,:,:) = carLM_pTT(:,nChan+1,:,:);
        
        carS_pTT(:,1+1:nChan+1,:,:) = carS_cTT(:,1:nChan,:,:);
        carS_pTT(:,1,:,:) = carS_pTT(:,nChan+1,:,:);
        
          %% information fiedelity
        chanCenter = ceil(nChan/2); % for shifted data, center is 4
        oriUnit = 2*pi/(nChan);
        oriRad = 0:oriUnit:2*pi;
        % oriRad = oriRad-oriRad(chanCenter);
        my_function = cos(abs(oriRad)-pi);
        
        
        % Cardinal vs. intercardinal
        for xCol = 1:nColCond
            for xAtt = 1:nAttCond
                xTT = pTT(:,:,xCol,xAtt);
                nnChan = size(xTT,2);
                c_center = ceil(nnChan/2);
                e = zeros(length(SN), nnChan);
                for i = 1:nnChan
                    e(:,i) = xTT(:,i) .* my_function(i);
                end
                perm_colSel{xPerm,1,1,xCol,xAtt} = mean(e,2);
                perm_colSel_SEM{xPerm,1,1,xCol,xAtt} =  std(perm_colSel{xPerm,1,1,xCol,xAtt},0,1) ...
                    ./ (sqrt(length(SN)-1));
            end
        end
        
        %SEM
        
        % L-M vs. S
        for xLMS = 1:nLMSCond
            for xAtt = 1:nAttCond
                car_xTT = car_pTT(:,:,xLMS,xAtt);
                nnChan = size(car_xTT,2);
                c_center = ceil(nnChan/2);
                e = zeros(length(SN), nnChan);
                for i = 1:nnChan
                    e(:,i) = car_xTT(:,i) .* my_function(i);
                end
                perm_colSel{xPerm,1,2,xLMS,xAtt} = mean(e,2);
                perm_colSel_SEM{xPerm,1,2,xLMS,xAtt} =  std(perm_colSel{xPerm,1,2,xLMS,xAtt},0,1) ...
                    ./ (sqrt(length(SN)-1));
            end
        end
        
        
        % L-M 1 vs. L-M 2
        for xLM = 1:nLMSCond
            for xAtt = 1:nAttCond
                carLM_xTT = carLM_pTT(:,:,xLM,xAtt);
                nnChan = size(carLM_xTT,2);
                c_center = ceil(nnChan/2);
                e = zeros(length(SN), nnChan);
                for i = 1:nnChan
                    e(:,i) = carLM_xTT(:,i) .* my_function(i);
                end
                perm_carLMS_colSel{xPerm,1,1,xLM,xAtt} = mean(e,2);
                perm_carLMS_colSel_SEM{xPerm,1,1,xLM,xAtt} =  std(perm_carLMS_colSel{xPerm,1,1,xLM,xAtt},0,1) ...
                    ./ (sqrt(length(SN)-1));
            end
        end
        
        
        % S 1 vs. S 2
        for xS = 1:nLMSCond
            for xAtt = 1:nAttCond
                carS_xTT = carS_pTT(:,:,xS,xAtt);
                nnChan = size(carS_xTT,2);
                c_center = ceil(nnChan/2);
                e = zeros(length(SN), nnChan);
                for i = 1:nnChan
                    e(:,i) = carS_xTT(:,i) .* my_function(i);
                end
                perm_carLMS_colSel{xPerm,1,2,xS,xAtt} = mean(e,2);
                perm_carLMS_colSel_SEM{xPerm,1,2,xS,xAtt} =  std(perm_carLMS_colSel{xPerm,1,2,xS,xAtt},0,1) ...
                    ./ (sqrt(length(SN)-1));
            end
        end
        
    end % perm
    
    
    %% load real data
     TT = cell(length(SN), 1); zTT = cell(length(SN), 1);
    for xSN = 1:length(SN)
        %         if isitLGN == 1
        %             if xSN == 2 || xSN == 6 || xSN ==9 || xSN == 10 || xSN == 1
        %                 ROIs{xROI} = 'LGN_hk3_p.05';
        %             else
        %                 ROIs{xROI} = 'LGN_hk2_p.05';
        %             end
        %         end
%         tunTxt = sprintf('_tuning_%s_2TR_used4TRs_shift', ROIs{xROI});
         tunTxt = sprintf('_tuning_%s2TR_used4TRs_neutralModel_shift', ROIs{xROI});
        fileName = sprintf('%s%s.txt', SN{xSN}, tunTxt);
        TT{xSN} = load(fullfile(baseDir, SN{xSN}, tunDir, fileName));
        
        % zscore across channel
        for rr = 1:nCond % cc: row, color x attention
            zTT{xSN}(rr,1:nChan) = zscore(TT{xSN}(rr, :));
        end
    end
    
    tTT = []; % zcTT(sub, channel, cond(color*attention)
    for xSN = 1:length(SN)
        for rr = 1:nCond
            tTT(xSN, 1:nChan, rr) = zTT{xSN}(rr,:);
        end
    end
    
    
    %% Tuning for each cond.
    % cTT(sub, channel, color, cond)
    cTT = [];
    cTT(:,1:nChan,1,1) = mean(tTT(:,:,carInx), 3); %cardinal, in
    cTT(:,1:nChan,2,1) = mean(tTT(:,:,intInx), 3); %intercardinal, in
    cTT(:,1:nChan,1,2) = mean(tTT(:,:,carInx+nColor), 3); %cardinal, out
    cTT(:,1:nChan,2,2) = mean(tTT(:,:,intInx+nColor), 3); %intercardinal, out
    car_cTT = [];
    car_cTT(:,1:nChan,1,1) = mean(tTT(:,:,LMInx), 3); %L-M color, in
    car_cTT(:,1:nChan,2,1) = mean(tTT(:,:,SInx), 3); %S color, in
    car_cTT(:,1:nChan,1,2) = mean(tTT(:,:,LMInx+nColor), 3); %L-M color, out
    car_cTT(:,1:nChan,2,2) = mean(tTT(:,:,SInx+nColor), 3); %S color, out
    carLM_cTT = [];
    carLM_cTT(:,1:nChan,1,1) = mean(tTT(:,:,LMInx(1)), 3); %L-M color1, in
    carLM_cTT(:,1:nChan,2,1) = mean(tTT(:,:,LMInx(2)), 3); %L-M color2, in
    carLM_cTT(:,1:nChan,1,2) = mean(tTT(:,:,LMInx(1)+nColor), 3); %L-M color1, out
    carLM_cTT(:,1:nChan,2,2) = mean(tTT(:,:,LMInx(2)+nColor), 3); %L-M color2, out
    carS_cTT = [];
    carS_cTT(:,1:nChan,1,1) = mean(tTT(:,:,SInx(1)), 3); %L-M color1, in
    carS_cTT(:,1:nChan,2,1) = mean(tTT(:,:,SInx(2)), 3); %L-M color2, in
    carS_cTT(:,1:nChan,1,2) = mean(tTT(:,:,SInx(1)+nColor), 3); %L-M color1, out
    carS_cTT(:,1:nChan,2,2) = mean(tTT(:,:,SInx(2)+nColor), 3); %L-M color2, out
    
    
    %SEM
    cSEM_CI = []; cSEM_LMS = [];
    cSEM_CI(1, 1:nChan,:,:) = std(cTT,0,1) ./ (sqrt(length(SN)-1));
    cSEM_LMS(1, 1:nChan,:,:) = std(car_cTT,0,1) ./ (sqrt(length(SN)-1));
    cSEM_carLM(1, 1:nChan,:,:) = std(carLM_cTT,0,1) ./ (sqrt(length(SN)-1));
    cSEM_carS(1, 1:nChan,:,:) = std(carS_cTT,0,1) ./ (sqrt(length(SN)-1));
    
    %% PLOT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % when shifted, tuning functions are centered on channel 4
    % channel 8 value is added before channel 1, so total of 9 channels exist.
    
    %% plot each subj's tuning
    
    % Add one more channel to the end to make it odd number
    % Add the last channel value as channel 9
    pTT = []; car_pTT = []; carLM_pTT = []; carS_pTT = [];
    
    pTT(:,1+1:nChan+1,:,:) = cTT(:,1:nChan,:,:);
    pTT(:,1,:,:) = pTT(:,nChan+1,:,:);
    
    car_pTT(:,1+1:nChan+1,:,:) = car_cTT(:,1:nChan,:,:);
    car_pTT(:,1,:,:) = car_pTT(:,nChan+1,:,:);
    
    carLM_pTT(:,1+1:nChan+1,:,:) = carLM_cTT(:,1:nChan,:,:);
    carLM_pTT(:,1,:,:) = carLM_pTT(:,nChan+1,:,:);
    
    carS_pTT(:,1+1:nChan+1,:,:) = carS_cTT(:,1:nChan,:,:);
    carS_pTT(:,1,:,:) = carS_pTT(:,nChan+1,:,:);
    
    
    
    colSel = cell(length(Methods), length(Hnames), 2,2);
    colSel_SEM = cell(length(Methods), length(Hnames), 2,2);
    carLMS_colSel = cell(length(Methods), length(Hnames), 2,2);
    carLMS_colSel_SEM = cell(length(Methods), length(Hnames), 2,2);
    
    
    % Cardinal vs. intercardinal
    for xCol = 1:nColCond
        for xAtt = 1:nAttCond
            xTT = pTT(:,:,xCol,xAtt);
            nnChan = size(xTT,2);
            c_center = ceil(nnChan/2);
            e = zeros(length(SN), nnChan);
            for i = 1:nnChan
                e(:,i) = xTT(:,i) .* my_function(i);
            end
            colSel{1,1,xCol,xAtt} = mean(e,2);
            colSel_SEM{1,1,xCol,xAtt} =  std(colSel{1,1,xCol,xAtt},0,1) ...
                ./ (sqrt(length(SN)-1));
        end
    end
    
    %SEM
    
    % L-M vs. S
    for xLMS = 1:nLMSCond
        for xAtt = 1:nAttCond
            car_xTT = car_pTT(:,:,xLMS,xAtt);
            nnChan = size(car_xTT,2);
            c_center = ceil(nnChan/2);
            e = zeros(length(SN), nnChan);
            for i = 1:nnChan
                e(:,i) = car_xTT(:,i) .* my_function(i);
            end
            colSel{1,2,xLMS,xAtt} = mean(e,2);
            colSel_SEM{1,2,xLMS,xAtt} =  std(colSel{1,2,xLMS,xAtt},0,1) ...
                ./ (sqrt(length(SN)-1));
        end
    end
    
    
    % L-M 1 vs. L-M 2
    for xLM = 1:nLMSCond
        for xAtt = 1:nAttCond
            carLM_xTT = carLM_pTT(:,:,xLM,xAtt);
            nnChan = size(carLM_xTT,2);
            c_center = ceil(nnChan/2);
            e = zeros(length(SN), nnChan);
            for i = 1:nnChan
                e(:,i) = carLM_xTT(:,i) .* my_function(i);
            end
            carLMS_colSel{1,1,xLM,xAtt} = mean(e,2);
            carLMS_colSel_SEM{1,1,xLM,xAtt} =  std(carLMS_colSel{1,1,xLM,xAtt},0,1) ...
                ./ (sqrt(length(SN)-1));
        end
    end
    
    
    % S 1 vs. S 2
    for xS = 1:nLMSCond
        for xAtt = 1:nAttCond
            carS_xTT = carS_pTT(:,:,xS,xAtt);
            nnChan = size(carS_xTT,2);
            c_center = ceil(nnChan/2);
            e = zeros(length(SN), nnChan);
            for i = 1:nnChan
                e(:,i) = carS_xTT(:,i) .* my_function(i);
            end
            carLMS_colSel{1,2,xS,xAtt} = mean(e,2);
            carLMS_colSel_SEM{1,2,xS,xAtt} =  std(carLMS_colSel{1,2,xS,xAtt},0,1) ...
                ./ (sqrt(length(SN)-1));
        end
    end
    
    
    %% reorganize all of the perm data for one ROI into one array  
    out_colSel = []; out_car_colSel = []; 
    out_carLM_colSel = []; out_carS_colSel =[];
    for xPerm=1:nPerm
        for xCol=1:2
            for xAtt=1:2
                %1:SN 2:perm No. 3:Col(1:car, 2:int) 4:Att(1:in, 2:out), 5:fidelity
                tmp_repSel = [];
                tmp_repSel(:,1) = 1:length(SN);
                tmp_repSel(:,2:4) = repmat([xPerm, xCol,xAtt], [length(SN),1]);
                tmp_repSel(:,5) = perm_colSel{xPerm,1,1,xCol,xAtt};
                out_colSel = [out_colSel; tmp_repSel]; %colsel of all conds for one ROI
                
                
                tmp_repSel = [];
                tmp_repSel(:,1) = 1:length(SN);
                tmp_repSel(:,2:4) = repmat([xPerm, xCol,xAtt], [length(SN),1]);
                tmp_repSel(:,5) = perm_colSel{xPerm,1,2,xCol,xAtt};
                out_car_colSel = [out_car_colSel; tmp_repSel]; %colsel of all conds for one ROI
                
                tmp_repSel = [];
                tmp_repSel(:,1) = 1:length(SN);
                tmp_repSel(:,2:4) = repmat([xPerm, xCol,xAtt], [length(SN),1]);
                tmp_repSel(:,5) = perm_carLMS_colSel{xPerm,1,1,xCol,xAtt};
                out_carLM_colSel = [out_carLM_colSel; tmp_repSel]; %colsel of all conds for one ROI
                
                tmp_repSel = [];
                tmp_repSel(:,1) = 1:length(SN);
                tmp_repSel(:,2:4) = repmat([xPerm, xCol,xAtt], [length(SN),1]);
                tmp_repSel(:,5) = perm_carLMS_colSel{xPerm,1,2,xCol,xAtt};
                out_carS_colSel = [out_carS_colSel; tmp_repSel]; %colsel of all conds for one ROI
            end
        end
    end
    
    %% add real data to the array as perm no. 1001
    for xPerm = 1001
        for xCol=1:2
            for xAtt=1:2
                %1:SN 2:perm No. 3:Col(1:car, 2:int) 4:Att(1:in, 2:out), 5:fidelity
                tmp_repSel = [];
                tmp_repSel(:,1) = 1:length(SN);
                tmp_repSel(:,2:4) = repmat([xPerm, xCol,xAtt], [length(SN),1]);
                tmp_repSel(:,5) = colSel{1,1,xCol,xAtt};
                out_colSel = [out_colSel; tmp_repSel]; %colsel of all conds for one ROI
                
                
                tmp_repSel = [];
                tmp_repSel(:,1) = 1:length(SN);
                tmp_repSel(:,2:4) = repmat([xPerm, xCol,xAtt], [length(SN),1]);
                tmp_repSel(:,5) = colSel{1,2,xCol,xAtt};
                out_car_colSel = [out_car_colSel; tmp_repSel]; %colsel of all conds for one ROI
                
                tmp_repSel = [];
                tmp_repSel(:,1) = 1:length(SN);
                tmp_repSel(:,2:4) = repmat([xPerm, xCol,xAtt], [length(SN),1]);
                tmp_repSel(:,5) = carLMS_colSel{1,1,xCol,xAtt};
                out_carLM_colSel = [out_carLM_colSel; tmp_repSel]; %colsel of all conds for one ROI
                
                tmp_repSel = [];
                tmp_repSel(:,1) = 1:length(SN);
                tmp_repSel(:,2:4) = repmat([xPerm, xCol,xAtt], [length(SN),1]);
                tmp_repSel(:,5) = carLMS_colSel{1,2,xCol,xAtt};
                out_carS_colSel = [out_carS_colSel; tmp_repSel]; %colsel of all conds for one ROI
            end
        end
    end
    
    writematrix(out_colSel, fullfile(statDir, sprintf('%s_perm_%s', ROIs{xROI}, statFileName)));
    writematrix(out_car_colSel, fullfile(statDir, sprintf('%s_perm_car_%s', ROIs{xROI}, statFileName)));
    writematrix(out_carLM_colSel, fullfile(statDir, sprintf('%s_perm_carLM_%s', ROIs{xROI}, statFileName)));
    writematrix(out_carS_colSel, fullfile(statDir, sprintf('%s_perm_carS_%s', ROIs{xROI}, statFileName)));

end


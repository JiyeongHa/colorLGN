% 200615 jyh
% plot tuning graphs of each subj & avg
% calculate selectivity index with three different methods
% : information fidelity, linear regression, & cosine similarity
% Statistical evaluation of color selectivity
% : one-sample t test with 0
% Repeated-measures ANOVA (color X attention)
clear all;


%% Experiment info. & directory

exp='v3';
% ROIs = {'V2_3dg_fmasked_q.001','V3_3dg_fmasked_q.001', 'V4v_fmasked_q.001' }; %, 'V2_3dg_fmasked_q.001','V3_3dg_fmasked_q.001', 'V4v_fmasked_q.001'};
% ROIs = {'V1_3dg_fmasked_q.001'}; %, 'V2_3dg_fmasked_q.001','V3_3dg_fmasked_q.001', 'V4v_fmasked_q.001'};
% ROIs = {'V1_3dg_fmasked_q.001'};

% ROIs = {'LGN_hk4_p.05', ,'V1_3dg_fmasked_q.001','V2_3dg_fmasked_q.001','V3_3dg_fmasked_q.001', 'V4v_fmasked_q.001'};
% ROIs = {'LGN_hk2_p.05'};
ROIs = {'V1_3dg_fmasked_q.001','V2_3dg_fmasked_q.001','V3_3dg_fmasked_q.001', 'V4v_fmasked_q.001'};

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
nPerm = 1000;


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
perm_tval = cell(length(ROIs), length(Methods), length(Hnames), 2,2);
perm_carLMS_tval = cell(length(ROIs), length(Methods), length(Hnames), 2,2);

perm_rmfit = zeros(nPerm, length(ROIs), length(Methods), length(Hnames), 3);  %main1,2,int1, %last dim: p,f, & sig value
perm_carLMS_rmfit = zeros(nPerm, length(ROIs), length(Methods), length(Hnames),3);
perm_fval = cell(length(ROIs), length(Methods), length(Hnames)); %main1, main2, int
perm_carLMS_fval = cell(length(ROIs), length(Methods), length(Hnames), 3);

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

%% write csv
csvDir = fullfile(baseDir, 'perm_fidelity_ptfvalues');
csvPrefix_pt = fullfile(csvDir, sprintf('%s_permed_fidelity_ptvalues_9chan_2TRlag_neutralModel.csv', exp));
csvPrefix_pf = fullfile(csvDir, sprintf('%s_permed_fidelity_pfvalues_9chan_2TRlag_neutralModel.csv', exp));

csv_pt = fopen(csvPrefix_pt, 'a+');
csv_pf = fopen(csvPrefix_pf, 'a+');

% write headers
fprintf(csv_pt, 'ROI, car-in, car-out, int-in, int-out, L-M-in, L-M-out, S-in, S-out\n');
fprintf(csv_pf, 'ROI, Fvalue, Pvalue\n');



%% load perm data

for xROI = 1:length(ROIs)
    

perm_colSel = cell(nPerm, length(Methods), length(Hnames), 2,2);
perm_carLMS_colSel = cell(nPerm, length(Methods), length(Hnames), 2,2);
perm_ttfit = zeros(nPerm, length(Methods), length(Hnames), 2, 2, 3); %last dim: p,t, & sig value
perm_carLMS_ttfit = zeros(nPerm, length(Methods), length(Hnames), 2, 2, 3); %last dim: p,t, & sig value
perm_ANOVA = cell(length(Methods), length(Hnames));

    fprintf('.....ROI: %s permutation start.....\n\n', ROIs{xROI})
    %% Subject info.
    if strcmp(ROIs{xROI}, 'LGN_hk') == 1
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
    else
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
        %         SN{end+1} = '12'; % LHB
        SN{end+1} = '13'; % JHY
        SN{end+1} = '14'; % KYJ
    end
    for xPerm = 1:nPerm
            clear TT zTT tTT cTT car_cTT carLM_cTT carS_cTT pTT car_pTT carLM_pTT carS_pTT; 
        if mod(xPerm,100) == 0; fprintf('ROI: %s testing no.%d\n', ROIs{xROI}, xPerm); end
        for xSN = 1:length(SN)
            %             if isitLGN == 1
            %                 if xSN == 2 || xSN == 6 || xSN ==9 || xSN == 10 || xSN == 1
            %                     ROIs{xROI} = 'LGN_hk3_p.05';
            %                 else
            %                     ROIs{xROI} = 'LGN_hk2_p.05';
            %                 end
            %             end
            % fileName = sprintf('perm%d_%s_tuning_%s_used4TRs_shift.txt', xPerm, SN{xSN}, ROIs{xROI});
            
            fileName = sprintf('perm%d_%s_tuning_%s2TR_used4TRs_neutralModel_shift.txt', xPerm, SN{xSN}, ROIs{xROI});
            
            TT{xSN} = load(fullfile(baseDir, SN{xSN}, permDir, [ROIs{xROI} '_2TRlag_neutralModel'], fileName));
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
        cSEM_CI = []; cSEM_LMS = []; cSEM_carLM = []; cSEM_carS = [];
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
        
        %% information fiedelity
        chanCenter = ceil(nChan/2); % for shifted data (pTT), center is 4
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
        
    end
    
    
    %% one-sample t-test with 0

    for xPerm = 1:nPerm
        fprintf('.....now permutation no.%d testing with 0.....\n', xPerm)
        for xMethods = 1:length(Methods)
            for xH = 1:length(Hnames) % Hypothesis to test: 1.CI, 2.L-Mvs.S
                for xCI = 1:nColCond %CI or L-MS
                    for xAtt = 1:nAttCond
                        [sig, p, ~, stat] = ttest(perm_colSel{xPerm,xMethods,xH,xCI,xAtt}); %sig ==1, diff from 0
                        perm_ttfit(xPerm, xMethods,xH,xCI,xAtt,1) = p;
                        perm_ttfit(xPerm, xMethods,xH,xCI,xAtt,2) = stat.tstat;
                        perm_ttfit(xPerm, xMethods,xH,xCI,xAtt,3) = sig;
                    end
                end
            end
        end
        
        for xMethods = 1:length(Methods)
            for xH = 1:length(Hnames) % Hypothesis to test: L-M 1 vs. 5, S 3 vs. 7
                for xCI = 1:nColCond %CI or L-MS
                    for xAtt = 1:nAttCond
                        [sig, p, ~, stat] = ttest(perm_carLMS_colSel{xPerm,xMethods,xH,xCI,xAtt}); %sig ==1, diff from 0
                        perm_carLMS_ttfit(xPerm, xMethods,xH,xCI,xAtt,1) = p;
                        perm_carLMS_ttfit(xPerm, xMethods,xH,xCI,xAtt,2) = stat.tstat;
                        perm_carLMS_ttfit(xPerm, xMethods,xH,xCI,xAtt,3) = sig;
                    end
                end
            end
        end
    end
    
    for xMethods = 1:length(Methods)
        for xH = 1:length(Hnames) % Hypothesis to test: L-M 1 vs. 5, S 3 vs. 7
            for xCI = 1:nColCond %CI or L-MS
                for xAtt = 1:nAttCond
                    perm_tval{xROI, xMethods, xH, xCI, xAtt} = perm_ttfit(:,xMethods,xH,xCI,xAtt,2);
                    perm_carLMS_tval{xROI, xMethods, xH, xCI, xAtt} = perm_carLMS_ttfit(:,xMethods,xH,xCI,xAtt,2);
                    
                end
            end
        end
    end
    
        %% Repeated-measures ANOVA
        for xPerm = 1:nPerm
            fprintf('.....now permutation no.%d testing ANOVA .....\n', xPerm)
            for xMethods = 1:length(Methods)
                for i = 1:length(Hnames) % 1 = all color, 2 = cardinal colors only
                    perm_ANOVA{xMethods,i}(:,1) = [perm_colSel{xPerm,xMethods,i,1,1}; perm_colSel{xPerm,xMethods,i,2,1}; perm_colSel{xPerm,xMethods,i,1,2}; perm_colSel{xPerm,xMethods,i,2,2}];
                    perm_ANOVA{xMethods,i}(:,2) = repmat([1:length(SN)]', [nColCond*nAttCond, 1]); % SN
                    perm_ANOVA{xMethods,i}(:,3) = repmat([ones(length(SN),1); ones(length(SN),1)*2], [nAttCond,1]); %factor 1: car vs. incar
                    perm_ANOVA{xMethods,i}(:,4) = [ones(length(SN)*2,1); ones(length(SN)*2,1)*2]; %factor2: in vs. out
    
                    xrmfit = rm_anova2(perm_ANOVA{xMethods,i}(:,1),perm_ANOVA{xMethods,i}(:,2),perm_ANOVA{xMethods,i}(:,3),perm_ANOVA{xMethods,i}(:,4),{factorNames{i,:}});
                    perm_rmfit(xPerm, xROI, xMethods, i, 1) = xrmfit{2,5}; % color main effect
                    perm_rmfit(xPerm, xROI, xMethods, i, 2) = xrmfit{3,5}; % attention main effect
                    perm_rmfit(xPerm, xROI, xMethods, i, 3) = xrmfit{4,5}; % interaction
    
%     
%                 t = table(perm_colSel{xPerm,xMethods,i,1,1}, perm_colSel{xPerm,xMethods,i,1,2} ,perm_colSel{xPerm,xMethods,i,2,1},perm_colSel{xPerm,xMethods,i,2,2},...
%                     'VariableNames', {'carin', 'carout', 'intin', 'intout'});
%                 WithinStructure = table([1 1 2 2]', [1 2 1 2]', 'VariableNames', {'Color', 'Attention'});
%                 WithinStructure.Color = categorical(WithinStructure.Color);
%                 WithinStructure.Attention = categorical(WithinStructure.Attention);
%                 %WithinStructure.Color_Attention = WithinStructure.Color .* WithinStructure.Attention;
%                 rm = fitrm(t, 'carin, carout, intin, intout ~ 1', 'WithinDesign', WithinStructure);
%                 rmtable = ranova(rm, 'WithinModel', 'Color*Attention');
%                 c = multcompare(rm, 'Attention', 'By', 'Color');
%                 a = multcompare(rm, 'Color', 'By', 'Attention');
%                 perm_multcomp_test(xPerm, xROI, xMethods, i, 1) = c.pValue(1); % when color1, attention effect
%                 perm_multcomp_test(xPerm, xROI, xMethods, i, 2) = c.pValue(3); % when color2, attention effect
%                 perm_multcomp_test(xPerm, xROI, xMethods, i, 3) = a.pValue(1); % when attention1, color effect
%                 perm_multcomp_test(xPerm, xROI, xMethods, i, 4) = a.pValue(3); % when attention2, color effect
    
    
                    perm_ANOVA{xMethods,i}(:,1) = [perm_carLMS_colSel{xPerm,xMethods,i,1,1}; perm_carLMS_colSel{xPerm,xMethods,i,2,1}; perm_carLMS_colSel{xPerm,xMethods,i,1,2}; perm_carLMS_colSel{xPerm,xMethods,i,2,2}];
                    perm_ANOVA{xMethods,i}(:,2) = repmat([1:length(SN)]', [nColCond*nAttCond, 1]); % SN
                    perm_ANOVA{1,i}(:,3) = repmat([ones(length(SN),1); ones(length(SN),1)*2], [nAttCond,1]); %factor 1: L-M 1 vs. L-M 5
                    perm_ANOVA{1,i}(:,4) = [ones(length(SN)*2,1); ones(length(SN)*2,1)*2]; %factor2: in vs. out
    
                    xrmfit = rm_anova2(perm_ANOVA{xMethods,i}(:,1),perm_ANOVA{xMethods,i}(:,2),perm_ANOVA{xMethods,i}(:,3),perm_ANOVA{xMethods,i}(:,4),{factorNames{i,:}});
                    perm_carLMS_rmfit(xPerm, xROI, xMethods, i, 1) = xrmfit{2,5}; % color main effect
                    perm_carLMS_rmfit(xPerm, xROI, xMethods, i, 2) = xrmfit{3,5}; % attention main effect
                    perm_carLMS_rmfit(xPerm, xROI, xMethods, i, 3) = xrmfit{4,5}; % interaction
%     
%     
%                 t = table(perm_carLMS_colSel{xPerm,xMethods,i,1,1}, perm_carLMS_colSel{xPerm,xMethods,i,1,2} ,perm_carLMS_colSel{xPerm,xMethods,i,2,1},perm_carLMS_colSel{xPerm,xMethods,i,2,2},...
%                     'VariableNames', {'carin', 'carout', 'intin', 'intout'});
%                 WithinStructure = table([1 1 2 2]', [1 2 1 2]', 'VariableNames', {'Color', 'Attention'});
%                 WithinStructure.Color = categorical(WithinStructure.Color);
%                 WithinStructure.Attention = categorical(WithinStructure.Attention);
%                 %WithinStructure.Color_Attention = WithinStructure.Color .* WithinStructure.Attention;
%                 rm = fitrm(t, 'carin, carout, intin, intout ~ 1', 'WithinDesign', WithinStructure);
%                 rmtable = ranova(rm, 'WithinModel', 'Color*Attention');
%                 c = multcompare(rm, 'Attention', 'By', 'Color');
%                 a = multcompare(rm, 'Color', 'By', 'Attention');
%                 perm_carLMS_multcomp_test(xPerm, xROI, xMethods, i, 1) = c.pValue(1); % when color1, attention effect
%                 perm_carLMS_multcomp_test(xPerm, xROI, xMethods, i, 2) = c.pValue(3); % when color2, attention effect
%                 perm_carLMS_multcomp_test(xPerm, xROI, xMethods, i, 3) = a.pValue(1); % when attention1, color effect
%                 perm_carLMS_multcomp_test(xPerm, xROI, xMethods, i, 4) = a.pValue(3); % when attention2, color effect
    
                end
            end
        end
    
    
    
    
end % ROI


%% load real data
for xROI = 1:length(ROIs)
            clear TT zTT tTT cTT car_cTT carLM_cTT carS_cTT pTT car_pTT carLM_pTT carS_pTT; 
    colSel = cell(nPerm, length(Methods), length(Hnames), 2,2);
    carLMS_colSel = cell(nPerm, length(Methods), length(Hnames), 2,2);
    ttfit = zeros(nPerm, length(Methods), length(Hnames), 2, 2, 3); %last dim: p,t, & sig value
    carLMS_ttfit = zeros(nPerm, length(Methods), length(Hnames), 2, 2, 3); %last dim: p,t, & sig value
    ANOVA = cell(length(Methods), length(Hnames));
    
    carLMS_colSel = cell(length(Methods), length(Hnames), 2,2);
    carLMS_colSel_SEM = cell(length(Methods), length(Hnames), 2,2);

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
    
    
    %% plot average tuning
    % average
    
    avg_CI = []; avg_LMS = []; avg_carLM = []; avg_carS = [];
    for xCI = 1:2
        for xAtt = 1:2
            avg_CI(1,1:nChan+1,xCI,xAtt) = mean(pTT(:,1:nChan+1,xCI,xAtt), 1); %cardinal vs. intercardinal
            avg_LMS(1,1:nChan+1,xCI,xAtt) = mean(car_pTT(:,1:nChan+1,xCI,xAtt), 1); %L-M vs. S
            avg_carLM(1,1:nChan+1,xCI,xAtt) = mean(carLM_pTT(:,1:nChan+1,xCI,xAtt), 1); %L-M 1, L-M 2
            avg_carS(1,1:nChan+1,xCI,xAtt) = mean(carS_pTT(:,1:nChan+1,xCI,xAtt), 1); %L-M 1, L-M 2
        end
    end
    
    pSEM_CI = []; pSEM_LMS = []; pSEM_carLM = []; pSEM_carS = [];
    %SEM
    pSEM_CI(1,1+1:nChan+1,:,:) = cSEM_CI(1,1:nChan,:,:);
    pSEM_CI(1,1,:,:) = pSEM_CI(1,nChan+1,:,:);
    pSEM_LMS(1,1+1:nChan+1,:,:) = cSEM_LMS(1,1:nChan,:,:);
    pSEM_LMS(1,1,:,:) = pSEM_LMS(1,nChan+1,:,:);
    pSEM_carLM(1,1+1:nChan+1,:,:) = cSEM_carLM(1,1:nChan,:,:);
    pSEM_carLM(1,1,:,:) = pSEM_carLM(1,nChan+1,:,:);
    pSEM_carS(1,1+1:nChan+1,:,:) = cSEM_carS(1,1:nChan,:,:);
    pSEM_carS(1,1,:,:) = pSEM_carS(1,nChan+1,:,:);
    
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
    
    %% one-sample t-test with 0
    sigtxt = {'not different', 'different'};
    
    fprintf('.....now testing real data with 0.....\n')
    for xMethods = 1:length(Methods)
        for xH = 1:length(Hnames) % Hypothesis to test: 1.CI, 2.L-Mvs.S
            for xCI = 1:nColCond %CI or L-MS
                for xAtt = 1:nAttCond
                    [sig, p, ~, stat] = ttest(colSel{xMethods,xH,xCI,xAtt}); %sig ==1, diff from 0
                    ttfit(xMethods,xH,xCI,xAtt,1) = p;
                    ttfit(xMethods,xH,xCI,xAtt,2) = stat.tstat;
                    ttfit(xMethods,xH,xCI,xAtt,3) = sig;
                end
            end
        end
    end
    
    
    for xMethods = 1:length(Methods)
        for xH = 1:length(Hnames) % Hypothesis to test: L-M 1 vs. 5, S 3 vs. 7
            for xCI = 1:nColCond %CI or L-MS
                for xAtt = 1:nAttCond
                    [sig, p, ~, stat] = ttest(carLMS_colSel{xMethods,xH,xCI,xAtt}); %sig ==1, diff from 0
                    carLMS_ttfit(xMethods,xH,xCI,xAtt,1) = p;
                    carLMS_ttfit(xMethods,xH,xCI,xAtt,2) = stat.tstat;
                    carLMS_ttfit(xMethods,xH,xCI,xAtt,3) = sig;
                end
            end
        end
    end
    
    for xMethods = 1:length(Methods)
        for xH = 1:length(Hnames) % Hypothesis to test: L-M 1 vs. 5, S 3 vs. 7
            for xCI = 1:nColCond %CI or L-MS
                for xAtt = 1:nAttCond
                    real_tval{xROI, xMethods, xH, xCI, xAtt} = ttfit(xMethods,xH,xCI,xAtt,2);
                    real_carLMS_tval{xROI, xMethods, xH, xCI, xAtt} = carLMS_ttfit(xMethods,xH,xCI,xAtt,2);
                    
                end
            end
        end
    end
    
        for xMethods = 1:length(Methods)
            for i = 1:length(Hnames) % 1 = all color, 2 = cardinal colors only
                ANOVA{xMethods,i}(:,1) = [colSel{xMethods,i,1,1}; colSel{xMethods,i,2,1}; colSel{xMethods,i,1,2}; colSel{xMethods,i,2,2}];
                ANOVA{xMethods,i}(:,2) = repmat([1:length(SN)]', [nColCond*nAttCond, 1]); % SN
                ANOVA{xMethods,i}(:,3) = repmat([ones(length(SN),1); ones(length(SN),1)*2], [nAttCond,1]); %factor 1: car vs. incar
                ANOVA{xMethods,i}(:,4) = [ones(length(SN)*2,1); ones(length(SN)*2,1)*2]; %factor2: in vs. out
    
                xrmfit = rm_anova2(ANOVA{xMethods,i}(:,1),ANOVA{xMethods,i}(:,2),ANOVA{xMethods,i}(:,3),ANOVA{xMethods,i}(:,4),{factorNames{i,:}});
                rmfit(xROI, xMethods, i, 1) = xrmfit{2,5}; % color main effect
                rmfit(xROI, xMethods, i, 2) = xrmfit{3,5}; % attention main effect
                rmfit(xROI, xMethods, i, 3) = xrmfit{4,5}; % interaction
%     
%                 % multiple comparison
%                 t = table(colSel{xMethods,i,1,1}, colSel{xMethods,i,1,2} ,colSel{xMethods,i,2,1},colSel{xMethods,i,2,2},...
%                     'VariableNames', {'carin', 'carout', 'intin', 'intout'});
%                 WithinStructure = table([1 1 2 2]', [1 2 1 2]', 'VariableNames', {'Color', 'Attention'});
%                 WithinStructure.Color = categorical(WithinStructure.Color);
%                 WithinStructure.Attention = categorical(WithinStructure.Attention);
%                 %WithinStructure.Color_Attention = WithinStructure.Color .* WithinStructure.Attention;
%                 rm = fitrm(t, 'carin, carout, intin, intout ~ 1', 'WithinDesign', WithinStructure);
%                 rmtable = ranova(rm, 'WithinModel', 'Color*Attention');
%                 c = multcompare(rm, 'Attention', 'By', 'Color');
%                 a = multcompare(rm, 'Color', 'By', 'Attention');
%                 multcomp_test(xROI, xMethods, i, 1) = c.pValue(1); % when color1, attention effect
%                 multcomp_test(xROI, xMethods, i, 2) = c.pValue(3); % when color2, attention effect
%                 multcomp_test(xROI, xMethods, i, 3) = a.pValue(1); % when attention1, color effect
%                 multcomp_test(xROI, xMethods, i, 4) = a.pValue(3); % when attention2, color effect
%     
%     
    
                ANOVA{xMethods,i}(:,1) = [carLMS_colSel{xMethods,i,1,1}; carLMS_colSel{xMethods,i,2,1}; carLMS_colSel{xMethods,i,1,2}; carLMS_colSel{xMethods,i,2,2}];
                ANOVA{xMethods,i}(:,2) = repmat([1:length(SN)]', [nColCond*nAttCond, 1]); % SN
                ANOVA{xMethods,i}(:,3) = repmat([ones(length(SN),1); ones(length(SN),1)*2], [nAttCond,1]); %factor 1: L-M 1 vs. L-M 5
                ANOVA{xMethods,i}(:,4) = [ones(length(SN)*2,1); ones(length(SN)*2,1)*2]; %factor2: in vs. out
    
                xrmfit = rm_anova2(ANOVA{xMethods,i}(:,1),ANOVA{xMethods,i}(:,2),ANOVA{xMethods,i}(:,3),ANOVA{xMethods,i}(:,4),{carLMS_factorNames{i,:}});
                carLMS_rmfit(xROI, xMethods, i, 1) = xrmfit{2,5}; % color main effect
                carLMS_rmfit(xROI, xMethods, i, 2) = xrmfit{3,5}; % attention main effect
                carLMS_rmfit(xROI, xMethods, i, 3) = xrmfit{4,5}; % interaction
    
%     
%     %                 % multiple comparison
%                 carLMS_t = table(colSel{xMethods,i,1,1}, colSel{xMethods,i,1,2} ,colSel{xMethods,i,2,1},colSel{xMethods,i,2,2},...
%                     'VariableNames', {'carin', 'carout', 'intin', 'intout'});
%                 carLMS_WithinStructure = table([1 1 2 2]', [1 2 1 2]', 'VariableNames', {'Color', 'Attention'});
%                 carLMS_WithinStructure.Color = categorical(carLMS_WithinStructure.Color);
%                 carLMS_WithinStructure.Attention = categorical(carLMS_WithinStructure.Attention);
%                 %WithinStructure.Color_Attention = WithinStructure.Color .* WithinStructure.Attention;
%                 carLMS_rm = fitrm(t, 'carin, carout, intin, intout ~ 1', 'WithinDesign', carLMS_WithinStructure);
%                 rmtable = ranova(carLMS_rm, 'WithinModel', 'Color*Attention');
%                 carLMS_c = multcompare(carLMS_rm, 'Attention', 'By', 'Color');
%                 carLMS_a = multcompare(carLMS_rm, 'Color', 'By', 'Attention');
%                 carLMS_multcomp_test(xROI, xMethods, i, 1) = carLMS_c.pValue(1); % when color1, attention effect
%                 carLMS_multcomp_test(xROI, xMethods, i, 2) = carLMS_c.pValue(3); % when color2, attention effect
%                 carLMS_multcomp_test(xROI, xMethods, i, 3) = carLMS_a.pValue(1); % when attention1, color effect
%                 carLMS_multcomp_test(xROI, xMethods, i, 4) = carLMS_a.pValue(3); % when attention2, color effect
%     
    
            end
        end
    
    
    
    
    
end % ROI

%% get P
FNL_tval = [];
FNL_fval = [];
FNL_pTval = [];
Effect = {'color', 'attention', 'interaction'};
ph_Effect = {'color1-attention effect', 'color2-attention effect', 'attention1-color effect', 'attention2-color effect'};

for xROI = 1:length(ROIs)
    for xMethods = 1:length(Methods)
        for xH = 1:length(Hnames) % all colors or car-colors
            if xH == 1, names = CInames; elseif xH == 2, names = LMSnames; end
            for xCI = 1:nColCond  % cardinal, intercardinal/ L-M, S
                for xAtt = 1:nAttCond % in or out
                    
                    xpermTval = perm_tval{xROI,xMethods,xH,xCI,xAtt}(1:nPerm);
                    xrealTval = real_tval{xROI,xMethods,xH,xCI,xAtt};
                    bigger_than_xrealTval = xpermTval(find(xrealTval <= xpermTval),:);
                    smaller_than_xrealTval = xpermTval(find(xrealTval >= xpermTval),:);
                    count_b_tval = (length(bigger_than_xrealTval)/nPerm);
                    count_s_tval = (length(smaller_than_xrealTval)/nPerm);
                    FNL_pval = count_b_tval;
                    sig = FNL_pval < 0.05;
                    fprintf('<%s> %s: %s_%s_%s is %s, t=%2.4f, p=%2.4f\n',...
                        ROIs{xROI}, Methods{xMethods}, Hnames{xH}, names{xCI}, Attnames{xAtt}, sigtxt{sig+1},...
                        xrealTval, min(count_b_tval, count_s_tval));
                    FNL_ptval(xROI,xMethods,xH,xCI,xAtt,1) =  FNL_pval;
                    FNL_ptval(xROI,xMethods,xH,xCI,xAtt,2) =  xrealTval;
                    
                    
                end
            end
        end
    end
    
    %write two rows for each ROI
    fprintf(csv_pt, '%s, %2.5f, %2.5f, %2.5f, %2.5f, %2.5f, %2.5f, %2.5f, %2.5f\n',...
        ROIs{xROI}, FNL_ptval(xROI,1,1,1,1,2),FNL_ptval(xROI,1,1,1,2,2),FNL_ptval(xROI,1,1,2,1,2),FNL_ptval(xROI,1,1,2,2,2),...
        FNL_ptval(xROI,1,2,1,1,2),FNL_ptval(xROI,1,2,1,2,2),FNL_ptval(xROI,1,2,2,1,2),FNL_ptval(xROI,1,2,2,2,2));
    fprintf(csv_pt, '%s, %2.5f, %2.5f, %2.5f, %2.5f, %2.5f, %2.5f, %2.5f, %2.5f\n',...
        ROIs{xROI}, FNL_ptval(xROI,1,1,1,1,1),FNL_ptval(xROI,1,1,1,2,1),FNL_ptval(xROI,1,1,2,1,1),FNL_ptval(xROI,1,1,2,2,1),...
        FNL_ptval(xROI,1,2,1,1,1),FNL_ptval(xROI,1,2,1,2,1),FNL_ptval(xROI,1,2,2,1,1),FNL_ptval(xROI,1,2,2,2,1));
    
    
end



for xROI = 1:length(ROIs)
    for xMethods = 1:length(Methods)
        for xH = 1:length(Hnames) % all colors or car-colors
            if xH == 1, names = CInames; elseif xH == 2, names = LMSnames; end
            for xEffect = 1:3 %main1, main2, interaction
                xpermFval = perm_rmfit(1:nPerm, xROI, xMethods, xH, xEffect);
                xrealFval = rmfit(xROI, xMethods, xH, xEffect);
                bigger_than_xrealFval = xpermFval(find(xrealFval <= xpermFval), :);
                FNL_ftest_pval = (length(bigger_than_xrealFval)/nPerm);
                sig = FNL_ftest_pval < 0.05;
                fprintf('<%s> %s: %s ANOVA %s is %s, F=%2.4f, p=%2.4f\n',...
                    ROIs{xROI}, Methods{xMethods}, Hnames{xH}, Effect{xEffect}, sigtxt{sig+1},...
                    xrealFval, FNL_ftest_pval);
                FNL_pfval(xROI, xMethods, xH, xEffect,1) = FNL_ftest_pval;
                FNL_pfval(xROI, xMethods, xH, xEffect,2) = xrealFval;
            end
        end
    end
    
    
    %write 6 rows for each ROI
    fprintf(csv_pf, '%s, %2.5f, %2.5f\n', ROIs{xROI}, FNL_pfval(xROI, 1,1,1,2),FNL_pfval(xROI, 1,1,1,1));
    fprintf(csv_pf, '%s, %2.5f, %2.5f\n', ROIs{xROI}, FNL_pfval(xROI, 1,1,2,2),FNL_pfval(xROI, 1,1,2,1));
    fprintf(csv_pf, '%s, %2.5f, %2.5f\n', ROIs{xROI}, FNL_pfval(xROI, 1,1,3,2),FNL_pfval(xROI, 1,1,3,1));
    fprintf(csv_pf, '%s, %2.5f, %2.5f\n', ROIs{xROI}, FNL_pfval(xROI, 1,2,1,2),FNL_pfval(xROI, 1,2,1,1));
    fprintf(csv_pf, '%s, %2.5f, %2.5f\n', ROIs{xROI}, FNL_pfval(xROI, 1,2,2,2),FNL_pfval(xROI, 1,2,2,1));
    fprintf(csv_pf, '%s, %2.5f, %2.5f\n', ROIs{xROI}, FNL_pfval(xROI, 1,2,3,2),FNL_pfval(xROI, 1,2,3,1));

    
    
    
end
%
% for xROI = 1:length(ROIs)
%     for xMethods = 1:length(Methods)
%         for xH = 1:length(Hnames) % all colors or car-colors
%             if xH == 1, names = CInames; elseif xH == 2, names = LMSnames; end
%             for xPh_Effect = 1:4  % cardinal, intercardinal/ L-M, S
%
%
%                     xperm_pTval = perm_multcomp_test(1:nPerm, xROI, xMethods, i, xPh_Effect);
%                     xreal_pTval = multcomp_test(xROI, xMethods, xH, xPh_Effect);
%                     bigger_than_xreal_pTval = xperm_pTval(find(xreal_pTval <= xperm_pTval),:);
%                     smaller_than_xreal_pTval = xperm_pTval(find(xreal_pTval >= xperm_pTval),:);
%                     count_b_tval = (length(bigger_than_xreal_pTval)/nPerm)*2;
%                     count_s_tval = (length(smaller_than_xreal_pTval)/nPerm)*2;
%                     sig = count_s_tval < 0.05;
%                     fprintf('<%s> %s: %s_%s is %s, t=%2.4f, p=%2.4f\n',...
%                         ROIs{xROI}, Methods{xMethods}, Hnames{xH}, ph_Effect{xPh_Effect}, sigtxt{sig+1},...
%                         xreal_pTval, count_s_tval);
%                     FNL_pTval(xROI,xMethods,xH,xPh_Effect,1) =  min(count_b_tval, count_s_tval);
%                     FNL_pTval(xROI,xMethods,xH,xPh_Effect,2) =  xreal_pTval;
%
%
%
%             end
%         end
%     end
% end
%
% %% histogram
% figure(1);
% for xMethods = 1:length(Methods)
%     for xH = 1:length(Hnames) % all colors or car-colors
%         if xH == 1, names = CInames; elseif xH == 2, names = LMSnames; end
%         for xCI = 1:nColCond  % cardinal, intercardinal/ L-M, S
%             for xAtt = 1:nAttCond % in or out
%                 subplot(length(Hnames),nColCond*nAttCond, (2*(xCI-1)+xAtt)+(xH-1)*4);
%                 histogram( perm_tval{xROI, xMethods, xH, xCI, xAtt}(1:nPerm,1), 30);
%                 hold on
%                 set(gca, 'box', 'off', 'TickDir', 'out', 'linewidth', 1.5);
%                 xlim([-8 8]); ylim([0 200]);
%                 if xH ==1 && xCI == 1 && xAtt == 1
%                     xline(real_tval{xROI,xMethods,xH,xCI,xAtt}, 'r', {'Actual t-value'});
%                 else
%                     xline(real_tval{xROI,xMethods,xH,xCI,xAtt}, 'r');
%
%                 end
%                 title(sprintf('%s %s', names{xCI}, Attnames{xAtt}));
%             end
%         end
%     end
% end
%
% %% histogram
% figure(2);
% for xMethods = 1:length(Methods)
%     for xH = 1:length(Hnames) % all colors or car-colors
%         for xPh_Effect = 1:4
%             subplot(length(Hnames),nColCond*nAttCond, xPh_Effect+4*(xH-1));
%             histogram(perm_multcomp_test(1:nPerm, xROI, xMethods, xH, xPh_Effect), 100);
%             hold on
%             set(gca, 'box', 'off', 'TickDir', 'out', 'linewidth', 1.5);
%             xlim([0 2]); ylim([0 50]);
%             if xH ==1 && xPh_Effect == 1
%                 xline(multcomp_test(xROI, xMethods, xH, xPh_Effect), 'r', {'Actual t-value'});
%             else
%                 xline(multcomp_test(xROI, xMethods, xH, xPh_Effect), 'r');
%             end
%         title(sprintf('%s', ph_Effect{xPh_Effect}));
%         end
%     end
% end
%



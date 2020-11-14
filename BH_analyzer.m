%% Behavioral data analyzer
addpath('/Users/auna/Script/');


%% Exp1: RSVP
% The number of targets presented is different in each trial. 
% (Which frame to present the target was determined by coinflip of .15
% probability.)

SN = {};

exp = 'RSVP';
bhFilePrefix = '_ColorTuningATTBlob_';

nColor = 8; %number of colors presented
nCond = 2; %spatial attention: in vs. out
nTrial = 32; % number of trials within a run
% nTarget = 9; % number of target frames within a trial
nCI = 2; % Cardinal vs. Intercardinal
nLMS = 2; % L-M vs. S
carInx = ([1:(nColor/2)]-1)*2+1; %1, 3, 5, 7
intInx = ([1:(nColor/2)])*2; %2,4,6,8
LMInx = [1,5];
SInx = [3,7];

nRun = 8;
nTR = 197;
blank1 = 2;
blank2 = 3;
sTR = 5;
fTR = 1;
tTR = sTR+fTR;
line_colors = {[1 0 0.6], [0.18 0.40 0.73]};

baseDir = sprintf('/Volumes/Duri/data/Color%s', exp);
bhDir = 'BH_data';


% plot parameters
lineWidth = 3;
lineColors = {[1 0.6 0.784],...
    [1 0 0.6],...
    [0.855 0.702 1],...
    [0.18 0.40 0.73],...
    [0 1 1],...
    [0 0.498 0],...
    [0.4 0.8 0],...
    [0.878 0.537 0.098]};
opacity = 0.3;

Hnames= {'All colors', 'Cardinal colors only'};
CInames = {'Cardinal', 'Int.cardinal'};
LMSnames = {'L-M', 'S'};
attNames = {'In', 'Out'};

attColors = {lineColors{2}, lineColors{4}};
CIColors = {lineColors{2}, lineColors{4}; lineColors{2}, lineColors{4}};
LMSColors = {lineColors{2}, lineColors{4}; lineColors{2}, lineColors{4}};

runNames = {'Run1', 'Run2', 'Run3', 'Run4', 'Run5', 'Run6', 'Run7', 'Run8'};
colorNames = {'Color1', 'Color2', 'Color3', 'Color4', 'Color5', 'Color6', 'Color7', 'Color8'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TotalAcc = [];
for xSN = 1:length(SN)
    for xRun = 1:nRun
        fileName = strcat(SN{xSN}{2}, bhFilePrefix, num2str(xRun), '_*.mat');
        fileList = dir(fullfile(baseDir, SN{xSN}{1}, bhDir, fileName));
        if length(fileList) > 1, error(sprintf('More than two data files for run %d exist!', xRun)), end
        
        load(fullfile(baseDir, SN{xSN}{1}, bhDir, fileList.name));
        targetAppear = saveInfo.RSVP; %frame length. nTrial X nTarget
        
        % The number of targets presented is different in each trial.
        % (Which frame to present the target was determined by coinflip of .15
        % probability.)
        % The frame in which target was presented is stored in saveInfo.RSVP(:,2)
        % and the number of frames presented for each quarter of a trial (reminder:
        % there were 4 parts in each trial, 4 quarter = 1 trial). So, we will
        % calculate how many frames were shown in each trial to track back targets
        % presented for each trial and store the trial info in saveInfo.RSVP(:,5).
        
        
        F4Q = saveInfo.Frame4Quarter;
        
        for i = 1:nTrial % how many frames for each trial (cumulated frames per trial)
            cF4T(i) = sum(F4Q(1:i*4));
        end
        
        % compare cumulated frames to the frame number for each target and store
        % the corresponding trial number in the fifth column.
        lastframe = 0;
        for i = 1:nTrial
            targetAppear(targetAppear(:,2)>lastframe & targetAppear(:,2)<=cF4T(i),5) = i;
            lastframe = cF4T(i);
        end

        for xTrial = 1:nTrial
            % trialInx = (xTrial-1)*nTarget +1:(xTrial*nTarget);
            trialAcc(xTrial,1) = mean(targetAppear(targetAppear(:,5) == xTrial, 4));
        end
        condMatrix = saveInfo.Condition; %1column: color, 2column:attention
        condAcc = [ones(nTrial,1)*xSN, ones(nTrial,1)*xRun, condMatrix trialAcc];
        TotalAcc = [TotalAcc; condAcc];
        
        for xColor = 1:nColor
            for xAtt = 1:nCond
                allSubCondAcc{xRun,xColor,xAtt,xSN} =...
                    condAcc(condAcc(:,1) == xColor & condAcc(:,2)==xAtt,3);
            end
        end
        
    end
end

if length(TotalAcc) ~= length(SN)*nRun*nTrial
    error('Some values are missing among runs, trials or subjects!')
end

%%TotalAcc
%%1SN 2Run 3Color 4Attention 5Accuracy

% Run acc.
% for each sub
runAcc = [];
for xSN =1:length(SN)
    for xRun = 1:nRun
        runAcc(xSN,xRun) = nanmean(TotalAcc(TotalAcc(:,1) == xSN & TotalAcc(:,2) == xRun, 5));
    end
    runSEM(1,1:xRun) = nanstd(runAcc,0,1) ./ sqrt(length(SN)-1);
end

% Color Acc.
colorAcc = [];
for xSN = 1:length(SN)
    for xColor = 1:nColor
        colorAcc(xSN,xColor) = nanmean(TotalAcc(TotalAcc(:,1) == xSN & TotalAcc(:,3) == xColor, 5));
    end
    colorSEM(1,1:xColor) = nanstd(colorAcc,0,1) ./ sqrt(length(SN)-1);
    
    
end

carColorAcc = mean(colorAcc(:,carInx),2);
carSEM = nanstd(reshape(colorAcc(:,carInx),[],1)) ./ sqrt(length(SN)-1);
intColorAcc = mean(colorAcc(:,intInx),2);
intSEM = nanstd(reshape(colorAcc(:,intInx),[],1)) ./ sqrt(length(SN)-1);
LMColorAcc = mean(colorAcc(:,LMInx),2);
LMSEM = nanstd(reshape(colorAcc(:,LMInx),[],1)) ./ sqrt(length(SN)-1);
SColorAcc = mean(colorAcc(:,SInx),2);
SSEM = nanstd(reshape(colorAcc(:,SInx),[],1)) ./ sqrt(length(SN)-1);

% Attention Acc.
attAcc = [];
for xSN = 1:length(SN)
    for xAtt = 1:nCond
        attAcc(xSN,xAtt) = nanmean(TotalAcc(TotalAcc(:,1) == xSN & TotalAcc(:,4) == xAtt, 5));
    end
    attSEM(1,1:xAtt) = nanstd(attAcc,0,1) ./ sqrt(length(SN)-1);
end


% Color vs. Attention (2 x 2) acc.
colorAttAcc = [];
for xSN = 1:length(SN)
    for xColor = 1:nColor
        for xAtt = 1:nCond
            colorAttAcc(xSN,xColor,xAtt) =...
                mean(TotalAcc(TotalAcc(:,1) == xSN &...
                TotalAcc(:,3) == xColor &...
                TotalAcc(:,4) == xAtt, 5));
        end
    end
end

CIAttAcc = [];
CIAttAcc(:,1,:) = mean(colorAttAcc(:,carInx,:),2); %bh acc. for Cardinal colors
CIAttAcc(:,2,:) = mean(colorAttAcc(:,intInx,:),2); %bh acc. for Intercardinal colors
LMSAttAcc = [];
LMSAttAcc(:,1,:) = mean(colorAttAcc(:,LMInx,:),2); %bh acc. for L-M colors
LMSAttAcc(:,2,:) = mean(colorAttAcc(:,SInx,:),2); %bh acc. for S colors
attSEM(1,1:xAtt) = nanstd(attAcc,0,1) ./ sqrt(length(SN)-1);


CIAttSEM = []; LMSAttSEM = [];
for xCI = 1:nCI
    for xAtt = 1:nCond
        CIAttSEM(1,xCI,xAtt) = nanstd([CIAttAcc(:,xCI,xAtt)],0,1) / sqrt(length(SN)-1);
        LMSAttSEM(1,xCI,xAtt) = nanstd([LMSAttAcc(:,xCI,xAtt)],0,1) / sqrt(length(SN)-1);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Figure 1: run
figure(1);
subplot(2, nRun, 1:nRun);
h = plot(nanmean(runAcc,1), 'color', 'k', 'linewidth', lineWidth);
hold on;
set(gca, 'box', 'off', 'XTickLabel', runNames);
xlim([1 8]); ylim([0 1]);
title(sprintf('%s color blob detection accuracy in each run', exp));
uSEM = []; dSEM = [];
uSEM = nanmean(runAcc,1)+runSEM(1,:);
dSEM = nanmean(runAcc,1)-runSEM(1,:);
yAxis = [uSEM fliplr(dSEM)];
xAxis = [1:length(uSEM) fliplr(1:length(dSEM))];
p = patch(xAxis, yAxis, get(h, 'Color'));
set(p, 'FaceAlpha', opacity, 'edgeColor', 'none');
xlim([1 8]); ylim([0 1]);

subplot(2, nColor, nRun+(1:nColor));
% Figure 2: color
b = bar(1:8, nanmean(colorAcc,1), 'EdgeColor', 'none');
hold on
title(sprintf('%s color blob detection accuracy in each color cond.', exp));
errorbar(1:8, mean(colorAcc,1), colorSEM, colorSEM, 'k', 'linestyle', 'none', 'LineWidth', 1);
set(gca, 'XTickLabel', colorNames, 'Box', 'off');
b.FaceColor = 'flat';
for xColor = 1:nColor
    b.CData(xColor,:) = lineColors{xColor};
end
ylim([0 1]);
hold off
set(gcf, 'renderer', 'painter');

% Figure 2: attention
figure(2);
b = bar(1:2, mean(attAcc,1), 'EdgeColor', 'none');
hold on
title(sprintf('%s color blob detection accuracy\n In vs. Out', exp));
errorbar(1:2, mean(attAcc,1), attSEM, attSEM, 'k', 'linestyle', 'none', 'LineWidth', 1);
set(gca, 'XTickLabel', attNames, 'Box', 'off');
b.FaceColor = 'flat';
for xAtt = 1:nCond
    b.CData(xAtt,:) = attColors{xAtt};
end
ylim([0 1]);
hold off
set(gcf, 'renderer', 'painter');

% Figure 3: color x attention
figure(3);

%Cardinal vs. Intercardinal
subplot(1,2,1);
color(1,:) = mean([CIAttAcc(:,1,1), CIAttAcc(:,1,2)],1); % e.g. cardinal: in,out
color(2,:) = mean([CIAttAcc(:,2,1), CIAttAcc(:,2,2)],1); % e.g. intercardinal: in,out
y = [color(1,:); color(2,:)];
sem(1,:) = [CIAttSEM(1,1,1) CIAttSEM(1,1,2)];
sem(2,:) = [CIAttSEM(1,2,1) CIAttSEM(1,2,2)];
err = [sem(1,:); sem(2,:)];

b = bar(y, 'EdgeColor', 'none');
hold on
for xAtt = 1:nCond
    b(xAtt).FaceColor = CIColors{1,xAtt};
end
set(gca, 'XTickLabel', CInames, 'Box', 'off');
title(sprintf('%s all colors x attention',exp));
[hh, icons, plots, txt] = legend(attNames, 'Location', 'Northeast');
nFactor1 = size(y,1);
nbars = size(y,2);
% Calculate the width for each bar group
barWidth = min(0.8, nbars/(nbars + 1.5));
% Set the position of each error bar in the center of the main bar
for i = 1:nbars
    x = (1:nFactor1) - barWidth/2 + (2*i-1)*barWidth/(2*nbars);
    errorbar(x, y(:,i), err(:,i), err(:,i), 'k', 'linestyle', 'none', 'LineWidth', 1);
end
ylim([0 1]);


%L-M vs. S
subplot(1,2,2);
color(1,:) = mean([LMSAttAcc(:,1,1), LMSAttAcc(:,1,2)],1); % e.g. cardinal: in,out
color(2,:) = mean([LMSAttAcc(:,2,1), LMSAttAcc(:,2,2)],1); % e.g. intercardinal: in,out
y = [color(1,:); color(2,:)];
sem(1,:) = [LMSAttSEM(1,1,1) LMSAttSEM(1,1,2)];
sem(2,:) = [LMSAttSEM(1,2,1) LMSAttSEM(1,2,2)];
err = [sem(1,:); sem(2,:)];

b = bar(y, 'EdgeColor', 'none');
hold on
for xAtt = 1:nCond
    b(xAtt).FaceColor = CIColors{1,xAtt};
end
set(gca, 'XTickLabel', LMSnames, 'Box', 'off');
nFactor1 = size(y,1);
nbars = size(y,2);
% Calculate the width for each bar group
barWidth = min(0.8, nbars/(nbars + 1.5));
% Set the position of each error bar in the center of the main bar
for i = 1:nbars
    x = (1:nFactor1) - barWidth/2 + (2*i-1)*barWidth/(2*nbars);
    errorbar(x, y(:,i), err(:,i), err(:,i), 'k', 'linestyle', 'none', 'LineWidth', 1);
end
ylim([0 1]);
title(sprintf('%s cardinal colors x attention',exp));
set(gcf, 'renderer', 'painter');


% RM ANOVA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('.\n.\n.\n.\n.\n.\n.\n.\n.\n.\n.\n.\n.\n');

% for run (one-way rm anova)
runRMtbl = array2table([[1:length(SN)]', runAcc]);
runRMtbl.Properties.VariableNames = ['SN' runNames];
within = table([1:nRun]', 'VariableNames', {'Run'});
runRMfit = fitrm(runRMtbl, 'Run1-Run8~1', 'WithinDesign', within);
rmANOVA.run = ranova(runRMfit);
if find(0.05 > rmANOVA.run.pValue), sigMsg = 'significant';
elseif find(rmANOVA.run.pValue > 0.05), sigMsg = 'not significant'; end
fprintf('************Run repeated-measures ANOVA results: %s***********\n', sigMsg);
rmANOVA.run

clear sigMsg;

% for colors (one-way rm anova)
colorRMtbl = array2table([[1:length(SN)]', colorAcc]);
colorRMtbl.Properties.VariableNames = ['SN' colorNames];
within = table([1:nColor]', 'VariableNames', {'Color'});
colorRMfit = fitrm(colorRMtbl, 'Color1-Color8~1', 'WithinDesign', within);
rmANOVA.color = ranova(colorRMfit);
if find(0.05 > rmANOVA.color.pValue), sigMsg = 'significant';
elseif find(rmANOVA.color.pValue > 0.05), sigMsg = 'not significant'; end
fprintf('************Color repeated-measures ANOVA results: %s***********\n', sigMsg);
rmANOVA.color
% posthoc.color = multcompare(colorRMfit, 'Color');
% posthoc.color = posthoc.color{find(posthoc.color.Color_1 < posthoc.color.Color_2),:};
fprintf('\n\n')


% for color x attention
% all colors
factorNames = cell(2,2);
factorNames = {'C vs. I', 'Attention'; 'L-M vs. S', 'Attention'};
CIanova = [];

CIanova(:,1) = reshape(CIAttAcc, [], 1);
CIanova(:,2) = repmat([1:length(SN)]', [nCI*nCond,1]);
CIanova(:,3) = repmat([ones(length(SN),1); ones(length(SN),1)*2], [nCI,1]); % factor1: CI
CIanova(:,4) = [ones(length(SN)*2,1); ones(length(SN)*2,1)*2]; % factor2
rmANOVA.CI = rm_anova2(CIanova(:,1), CIanova(:,2), CIanova(:,3), CIanova(:,4), {factorNames{1,:}});
fprintf('************ %s X attention ANOVA results ************\n', Hnames{1})
rmANOVA.CI
fprintf('\n')
LMSanova = [];

LMSanova(:,1) = reshape(LMSAttAcc, [], 1);
LMSanova(:,2) = repmat([1:length(SN)]', [nLMS*nCond,1]);
LMSanova(:,3) = repmat([ones(length(SN),1); ones(length(SN),1)*2], [nLMS,1]); % factor1: LMS
LMSanova(:,4) = [ones(length(SN)*2,1); ones(length(SN)*2,1)*2]; % factor2
rmANOVA.LMS = rm_anova2(LMSanova(:,1), LMSanova(:,2), LMSanova(:,3), LMSanova(:,4), {factorNames{2,:}});
fprintf('************ %s X attention ANOVA results ************\n', Hnames{2})
rmANOVA.LMS

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save .csv files for JAMOVI...

csvPrefixRun = fullfile(csvDir, sprintf('%s_datRunAcc.csv', exp));
csvPrefixColor = fullfile(csvDir, sprintf('%s_datColorAcc.csv', exp));
csvPrefixCI = fullfile(csvDir, sprintf('%s_datCIAttAcc.csv', exp));
csvPrefixLMS = fullfile(csvDir, sprintf('%s_datLMSAttAcc.csv', exp));

csvRun = fopen(csvPrefixRun, 'a+');
csvColor = fopen(csvPrefixColor, 'a+');
csvCI = fopen(csvPrefixCI, 'a+');
csvLMS = fopen(csvPrefixLMS, 'a+');

% write headers
fprintf(csvRun, 'Run1, Run2, Run3, Run4, Run5, Run6, Run7, Run8\n');
fprintf(csvColor, 'Color1, Color2, Color3, Color4, Color5, Color6, Color7, Color8\n');
fprintf(csvCI, 'C-In, C-Out, I-In, I-Out\n');
fprintf(csvLMS, 'LM-In, S-Out, LM-In, S-Out\n');

% write accuracies
for xSN = 1:length(SN)

fprintf(csvRun, '%2.4f, %2.4f, %2.4f, %2.4f, %2.4f, %2.4f, %2.4f, %2.4f\n',...
    runAcc(xSN,1), runAcc(xSN,2), runAcc(xSN,3), runAcc(xSN,4),...
    runAcc(xSN,5), runAcc(xSN,6), runAcc(xSN,7), runAcc(xSN,8));

fprintf(csvColor, '%2.4f, %2.4f, %2.4f, %2.4f, %2.4f, %2.4f, %2.4f, %2.4f\n',...
    colorAcc(xSN,1), colorAcc(xSN,2), colorAcc(xSN,3), colorAcc(xSN,4),...
    colorAcc(xSN,5), colorAcc(xSN,6), colorAcc(xSN,7), colorAcc(xSN,8));

fprintf(csvCI, '%2.4f, %2.4f, %2.4f, %2.4f\n',...
    CIAttAcc(xSN,1,1), CIAttAcc(xSN,1,2), CIAttAcc(xSN,2,1), CIAttAcc(xSN,2,2)); 

fprintf(csvLMS, '%2.4f, %2.4f, %2.4f, %2.4f\n',...
    LMSAttAcc(xSN,1,1), LMSAttAcc(xSN,1,2), LMSAttAcc(xSN,2,1), LMSAttAcc(xSN,2,2)); 

end

%% Exp2: Ecc
clear all;

SN={};
% SN{end+1} = {'01','KB'};
% SN{end+1} = {'02','JYM'};
% % SN{end+1} = {'03','KIS'};
% % SN{end+1} = {'04','SHY'};
SN{end+1} = {'05','JYS'};
SN{end+1} = {'06','PJH'};
SN{end+1} = {'07','PSY'};
SN{end+1} = {'08','JYJ'};
SN{end+1} = {'09','KDH'};
SN{end+1} = {'10','LYR'};
% SN{end+1} = {'11','SHY'};
SN{end+1} = {'12','KIS'};
SN{end+1} = {'13','SSK'};
SN{end+1} = {'14','HSH'};
SN{end+1} = {'16','BKD'};
SN{end+1} = {'15','HJY'};

exp = 'Ecc';
bhFilePrefix = '_ColorTuningATTEcc_';

nColor = 8; %number of colors presented
nCond = 2; %spatial attention: in vs. out
nTrial = 32; % number of trials within a run
nTarget = 9; % number of target frames within a trial
nCI = 2; % Cardinal vs. Intercardinal
nLMS = 2; % L-M vs. S
carInx = ([1:(nColor/2)]-1)*2+1; %1, 3, 5, 7
intInx = ([1:(nColor/2)])*2; %2,4,6,8
LMInx = [1,5];
SInx = [3,7];

nRun = 8;
nTR = 197;
blank1 = 2;
blank2 = 3;
sTR = 5;
fTR = 1;
tTR = sTR+fTR;
line_colors = {[1 0 0.6], [0.18 0.40 0.73]};

baseDir = sprintf('/group_hpc/WMShimLab2/PSY_Color/Color%s', exp);
bhDir = 'BH_data';
csvDir = '/home/jiyeongha/';

% plot parameters
lineWidth = 3;
lineColors = {[1 0.6 0.784],...
    [1 0 0.6],...
    [0.855 0.702 1],...
    [0.18 0.40 0.73],...
    [0 1 1],...
    [0 0.498 0],...
    [0.4 0.8 0],...
    [0.878 0.537 0.098]};
opacity = 0.3;

Hnames= {'All colors', 'Cardinal colors only'};
CInames = {'Cardinal', 'Int.cardinal'};
LMSnames = {'L-M', 'S'};
attNames = {'In', 'Out'};

attColors = {lineColors{2}, lineColors{4}};
CIColors = {lineColors{2}, lineColors{4}; lineColors{2}, lineColors{4}};
LMSColors = {lineColors{2}, lineColors{4}; lineColors{2}, lineColors{4}};

runNames = {'Run1', 'Run2', 'Run3', 'Run4', 'Run5', 'Run6', 'Run7', 'Run8'};
colorNames = {'Color1', 'Color2', 'Color3', 'Color4', 'Color5', 'Color6', 'Color7', 'Color8'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

TotalAcc = [];
for xSN = 1:length(SN)-1
    for xRun = 1:nRun
        fileName = strcat(SN{xSN}{2}, bhFilePrefix, num2str(xRun), '*.mat');
        if strcmp(SN{xSN}{2},'KIS') || strcmp(SN{xSN}{2},'PSY')
            fileName = strcat(SN{xSN}{1}, SN{xSN}{2}, bhFilePrefix, num2str(xRun), '*.mat');
        end
        fileList = dir(fullfile(baseDir, SN{xSN}{1}, bhDir, fileName));
        if length(fileList) > 1, error(sprintf('More than two data files for run %d exist!', xRun)), end
        
        load(fullfile(baseDir, SN{xSN}{1}, bhDir, fileList.name));
        targetAppear = saveInfo.RSVP; %frame length. nTrial X nTarget
        
        for xTrial = 1:nTrial
            trialInx = (xTrial-1)*nTarget +1:(xTrial*nTarget);
            trialAcc(xTrial,1) = mean(targetAppear(trialInx, 5));
        end
        condMatrix = saveInfo.Condition; %1column: color, 2column:attention
        condAcc = [ones(nTrial,1)*xSN, ones(nTrial,1)*xRun, condMatrix trialAcc];
        TotalAcc = [TotalAcc; condAcc];
        
        for xColor = 1:nColor
            for xAtt = 1:nCond
                allSubCondAcc{xRun,xColor,xAtt,xSN} =...
                    condAcc(condAcc(:,1) == xColor & condAcc(:,2)==xAtt,3);
            end
        end
        
    end
end

HJY_run = [1,5,6,7,8];
for xSN = 11
    for xRun = 1:length(HJY_run)
        fileName = strcat(SN{xSN}{2}, bhFilePrefix, num2str(HJY_run(xRun)), '*.mat');
        
        fileList = dir(fullfile(baseDir, SN{xSN}{1}, bhDir, fileName));
        if length(fileList) > 1, error(sprintf('More than two data files for run %d exist!', xRun)), end
        
        load(fullfile(baseDir, SN{xSN}{1}, bhDir, fileList.name));
        targetAppear = saveInfo.RSVP; %frame length. nTrial X nTarget
        
        for xTrial = 1:nTrial
            trialInx = (xTrial-1)*nTarget +1:(xTrial*nTarget);
            trialAcc(xTrial,1) = mean(targetAppear(trialInx, 5));
        end
        condMatrix = saveInfo.Condition; %1column: color, 2column:attention
        condAcc = [ones(nTrial,1)*xSN, ones(nTrial,1)*HJY_run(xRun), condMatrix trialAcc];
        TotalAcc = [TotalAcc; condAcc];
        
    end
end

if length(TotalAcc) ~= length(SN)*nRun*nTrial
    fprintf('Some values are missing among runs, trials or subjects!\n')
end

%%TotalAcc
%%1SN 2Run 3Color 4Attention 5Accuracy

% Run acc.
% for each sub
runAcc = [];
for xSN =1:length(SN)
    for xRun = 1:nRun
        runAcc(xSN,xRun) = nanmean(TotalAcc(TotalAcc(:,1) == xSN & TotalAcc(:,2) == xRun, 5));
    end
    runSEM(1,1:xRun) = nanstd(runAcc,0,1) ./ sqrt(length(SN)-1);
end

% Color Acc.
colorAcc = [];
for xSN = 1:length(SN)
    for xColor = 1:nColor
        colorAcc(xSN,xColor) = nanmean(TotalAcc(TotalAcc(:,1) == xSN & TotalAcc(:,3) == xColor, 5));
    end
    colorSEM(1,1:xColor) = nanstd(colorAcc,0,1) ./ sqrt(length(SN)-1);
    
    
end

carColorAcc = mean(colorAcc(:,carInx),2);
carSEM = nanstd(reshape(colorAcc(:,carInx),[],1)) ./ sqrt(length(SN)-1);
intColorAcc = mean(colorAcc(:,intInx),2);
intSEM = nanstd(reshape(colorAcc(:,intInx),[],1)) ./ sqrt(length(SN)-1);
LMColorAcc = mean(colorAcc(:,LMInx),2);
LMSEM = nanstd(reshape(colorAcc(:,LMInx),[],1)) ./ sqrt(length(SN)-1);
SColorAcc = mean(colorAcc(:,SInx),2);
SSEM = nanstd(reshape(colorAcc(:,SInx),[],1)) ./ sqrt(length(SN)-1);

% Attention Acc.
attAcc = [];
for xSN = 1:length(SN)
    for xAtt = 1:nCond
        attAcc(xSN,xAtt) = nanmean(TotalAcc(TotalAcc(:,1) == xSN & TotalAcc(:,4) == xAtt, 5));
    end
    attSEM(1,1:xAtt) = nanstd(attAcc,0,1) ./ sqrt(length(SN)-1);
end

% Color vs. Attention (2 x 2) acc.
colorAttAcc = [];
for xSN = 1:length(SN)
    for xColor = 1:nColor
        for xAtt = 1:nCond
            colorAttAcc(xSN,xColor,xAtt) =...
                mean(TotalAcc(TotalAcc(:,1) == xSN &...
                TotalAcc(:,3) == xColor &...
                TotalAcc(:,4) == xAtt, 5));
        end
    end
end

CIAttAcc = [];
CIAttAcc(:,1,:) = mean(colorAttAcc(:,carInx,:),2); %bh acc. for Cardinal colors
CIAttAcc(:,2,:) = mean(colorAttAcc(:,intInx,:),2); %bh acc. for Intercardinal colors
LMSAttAcc = [];
LMSAttAcc(:,1,:) = mean(colorAttAcc(:,LMInx,:),2); %bh acc. for L-M colors
LMSAttAcc(:,2,:) = mean(colorAttAcc(:,SInx,:),2); %bh acc. for S colors


CIAttSEM = []; LMSAttSEM = [];
for xCI = 1:nCI
    for xAtt = 1:nCond
        CIAttSEM(1,xCI,xAtt) = nanstd([CIAttAcc(:,xCI,xAtt)],0,1) / sqrt(length(SN)-1);
        LMSAttSEM(1,xCI,xAtt) = nanstd([LMSAttAcc(:,xCI,xAtt)],0,1) / sqrt(length(SN)-1);
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Figure 1: run
figure(1);
subplot(2, nRun, 1:nRun);
h = plot(nanmean(runAcc,1), 'color', 'k', 'linewidth', lineWidth);
hold on;
set(gca, 'box', 'off', 'XTickLabel', runNames);
xlim([1 8]); ylim([0 1]);
title(sprintf('%s color blob detection accuracy in each run', exp));
uSEM = []; dSEM = [];
uSEM = nanmean(runAcc,1)+runSEM(1,:);
dSEM = nanmean(runAcc,1)-runSEM(1,:);
yAxis = [uSEM fliplr(dSEM)];
xAxis = [1:length(uSEM) fliplr(1:length(dSEM))];
p = patch(xAxis, yAxis, get(h, 'Color'));
set(p, 'FaceAlpha', opacity, 'edgeColor', 'none');
xlim([1 8]); ylim([0 1]);

subplot(2, nColor, nRun+(1:nColor));
% Figure 2: color
b = bar(1:8, nanmean(colorAcc,1), 'EdgeColor', 'none');
hold on
title(sprintf('%s color blob detection accuracy in each color cond.', exp));
errorbar(1:8, mean(colorAcc,1), colorSEM, colorSEM, 'k', 'linestyle', 'none', 'LineWidth', 1);
set(gca, 'XTickLabel', colorNames, 'Box', 'off');
b.FaceColor = 'flat';
for xColor = 1:nColor
    b.CData(xColor,:) = lineColors{xColor};
end
ylim([0 1]);
hold off
set(gcf, 'renderer', 'painter');

% Figure 2: attention
figure(2);
b = bar(1:2, mean(attAcc,1), 'EdgeColor', 'none');
hold on
title(sprintf('%s color blob detection accuracy\n In vs. Out', exp));
errorbar(1:2, mean(attAcc,1), attSEM, attSEM, 'k', 'linestyle', 'none', 'LineWidth', 1);
set(gca, 'XTickLabel', attNames, 'Box', 'off');
b.FaceColor = 'flat';
for xAtt = 1:nCond
    b.CData(xAtt,:) = attColors{xAtt};
end
ylim([0 1]);
hold off
set(gcf, 'renderer', 'painter');

% Figure 3: color x attention
figure(3);

%Cardinal vs. Intercardinal
subplot(1,2,1);
color(1,:) = mean([CIAttAcc(:,1,1), CIAttAcc(:,1,2)],1); % e.g. cardinal: in,out
color(2,:) = mean([CIAttAcc(:,2,1), CIAttAcc(:,2,2)],1); % e.g. intercardinal: in,out
y = [color(1,:); color(2,:)];
sem(1,:) = [CIAttSEM(1,1,1) CIAttSEM(1,1,2)];
sem(2,:) = [CIAttSEM(1,2,1) CIAttSEM(1,2,2)];
err = [sem(1,:); sem(2,:)];

b = bar(y, 'EdgeColor', 'none');
hold on
for xAtt = 1:nCond
    b(xAtt).FaceColor = CIColors{1,xAtt};
end
set(gca, 'XTickLabel', CInames, 'Box', 'off');
title(sprintf('%s all colors x attention',exp));
[hh, icons, plots, txt] = legend(attNames, 'Location', 'Northeast');
nFactor1 = size(y,1);
nbars = size(y,2);
% Calculate the width for each bar group
barWidth = min(0.8, nbars/(nbars + 1.5));
% Set the position of each error bar in the center of the main bar
for i = 1:nbars
    x = (1:nFactor1) - barWidth/2 + (2*i-1)*barWidth/(2*nbars);
    errorbar(x, y(:,i), err(:,i), err(:,i), 'k', 'linestyle', 'none', 'LineWidth', 1);
end
ylim([0 1]);

%L-M vs. S
subplot(1,2,2);
color(1,:) = mean([LMSAttAcc(:,1,1), LMSAttAcc(:,1,2)],1); % e.g. cardinal: in,out
color(2,:) = mean([LMSAttAcc(:,2,1), LMSAttAcc(:,2,2)],1); % e.g. intercardinal: in,out
y = [color(1,:); color(2,:)];
sem(1,:) = [LMSAttSEM(1,1,1) LMSAttSEM(1,1,2)];
sem(2,:) = [LMSAttSEM(1,2,1) LMSAttSEM(1,2,2)];
err = [sem(1,:); sem(2,:)];

b = bar(y, 'EdgeColor', 'none');
hold on
for xAtt = 1:nCond
    b(xAtt).FaceColor = CIColors{1,xAtt};
end
set(gca, 'XTickLabel', LMSnames, 'Box', 'off');
nFactor1 = size(y,1);
nbars = size(y,2);
% Calculate the width for each bar group
barWidth = min(0.8, nbars/(nbars + 1.5));
% Set the position of each error bar in the center of the main bar
for i = 1:nbars
    x = (1:nFactor1) - barWidth/2 + (2*i-1)*barWidth/(2*nbars);
    errorbar(x, y(:,i), err(:,i), err(:,i), 'k', 'linestyle', 'none', 'LineWidth', 1);
end
ylim([0 1]);
title(sprintf('%s cardinal colors x attention',exp));
set(gcf, 'renderer', 'painter');


% RM ANOVA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('.\n.\n.\n.\n.\n.\n.\n.\n.\n.\n.\n.\n.\n');

% for run (one-way rm anova)
runRMtbl = array2table([[1:length(SN)]', runAcc]);
runRMtbl.Properties.VariableNames = ['SN' runNames];
within = table([1:nRun]', 'VariableNames', {'Run'});
runRMfit = fitrm(runRMtbl, 'Run1-Run8~1', 'WithinDesign', within);
rmANOVA.run = ranova(runRMfit);
if find(0.05 > rmANOVA.run.pValue), sigMsg = 'significant';
elseif find(rmANOVA.run.pValue > 0.05), sigMsg = 'not significant'; end
fprintf('************Run repeated-measures ANOVA results: %s***********\n', sigMsg);
rmANOVA.run

clear sigMsg;

% for colors (one-way rm anova)
colorRMtbl = array2table([[1:length(SN)]', colorAcc]);
colorRMtbl.Properties.VariableNames = ['SN' colorNames];
within = table([1:nColor]', 'VariableNames', {'Color'});
colorRMfit = fitrm(colorRMtbl, 'Color1-Color8~1', 'WithinDesign', within);
rmANOVA.color = ranova(colorRMfit);
if find(0.05 > rmANOVA.color.pValue), sigMsg = 'significant';
elseif find(rmANOVA.color.pValue > 0.05), sigMsg = 'not significant'; end
fprintf('************Color repeated-measures ANOVA results: %s***********\n', sigMsg);
rmANOVA.color
% posthoc.color = multcompare(colorRMfit, 'Color');
% posthoc.color = posthoc.color{find(posthoc.color.Color_1 < posthoc.color.Color_2),:};
fprintf('\n\n')


% for color x attention
% all colors
factorNames = cell(2,2);
factorNames = {'C vs. I', 'Attention'; 'L-M vs. S', 'Attention'};

CIanova(:,1) = reshape(CIAttAcc, [], 1);
CIanova(:,2) = repmat([1:length(SN)]', [nCI*nCond,1]);
CIanova(:,3) = repmat([ones(length(SN),1); ones(length(SN),1)*2], [nCI,1]); % factor1: CI
CIanova(:,4) = [ones(length(SN)*2,1); ones(length(SN)*2,1)*2]; % factor2
rmANOVA.CI = rm_anova2(CIanova(:,1), CIanova(:,2), CIanova(:,3), CIanova(:,4), {factorNames{1,:}});
fprintf('************ %s X attention ANOVA results ************\n', Hnames{1})
rmANOVA.CI
fprintf('\n')

LMSanova(:,1) = reshape(LMSAttAcc, [], 1);
LMSanova(:,2) = repmat([1:length(SN)]', [nLMS*nCond,1]);
LMSanova(:,3) = repmat([ones(length(SN),1); ones(length(SN),1)*2], [nLMS,1]); % factor1: LMS
LMSanova(:,4) = [ones(length(SN)*2,1); ones(length(SN)*2,1)*2]; % factor2
rmANOVA.LMS = rm_anova2(LMSanova(:,1), LMSanova(:,2), LMSanova(:,3), LMSanova(:,4), {factorNames{2,:}});
fprintf('************ %s X attention ANOVA results ************\n', Hnames{2})
rmANOVA.LMS

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save .csv files for JAMOVI...

csvPrefixRun = fullfile(csvDir, sprintf('%s_datRunAcc.csv', exp));
csvPrefixColor = fullfile(csvDir, sprintf('%s_datColorAcc.csv', exp));
csvPrefixCI = fullfile(csvDir, sprintf('%s_datCIAttAcc.csv', exp));
csvPrefixLMS = fullfile(csvDir, sprintf('%s_datLMSAttAcc.csv', exp));

csvRun = fopen(csvPrefixRun, 'a+');
csvColor = fopen(csvPrefixColor, 'a+');
csvCI = fopen(csvPrefixCI, 'a+');
csvLMS = fopen(csvPrefixLMS, 'a+');

% write headers
fprintf(csvRun, 'Run1, Run2, Run3, Run4, Run5, Run6, Run7, Run8\n');
fprintf(csvColor, 'Color1, Color2, Color3, Color4, Color5, Color6, Color7, Color8\n');
fprintf(csvCI, 'C-In, C-Out, I-In, I-Out\n');
fprintf(csvLMS, 'LM-In, S-Out, LM-In, S-Out\n');

% write accuracies
for xSN = 1:length(SN)

fprintf(csvRun, '%2.4f, %2.4f, %2.4f, %2.4f, %2.4f, %2.4f, %2.4f, %2.4f\n',...
    runAcc(xSN,1), runAcc(xSN,2), runAcc(xSN,3), runAcc(xSN,4),...
    runAcc(xSN,5), runAcc(xSN,6), runAcc(xSN,7), runAcc(xSN,8));

fprintf(csvColor, '%2.4f, %2.4f, %2.4f, %2.4f, %2.4f, %2.4f, %2.4f, %2.4f\n',...
    colorAcc(xSN,1), colorAcc(xSN,2), colorAcc(xSN,3), colorAcc(xSN,4),...
    colorAcc(xSN,5), colorAcc(xSN,6), colorAcc(xSN,7), colorAcc(xSN,8));

fprintf(csvCI, '%2.4f, %2.4f, %2.4f, %2.4f\n',...
    CIAttAcc(xSN,1,1), CIAttAcc(xSN,1,2), CIAttAcc(xSN,2,1), CIAttAcc(xSN,2,2)); 

fprintf(csvLMS, '%2.4f, %2.4f, %2.4f, %2.4f\n',...
    LMSAttAcc(xSN,1,1), LMSAttAcc(xSN,1,2), LMSAttAcc(xSN,2,1), LMSAttAcc(xSN,2,2)); 

end


%% Exp3: v3
clear all;

% subject info.
SN={};
SN{end+1} = {'01','PSY'};
SN{end+1} = {'03','KIS'};
SN{end+1} = {'04','JBH'};
SN{end+1} = {'05','HJH'};
SN{end+1} = {'06','KB'};
SN{end+1} = {'07','SHY'};
SN{end+1} = {'08','LSY'};
SN{end+1} = {'09','HJH'};
SN{end+1} = {'10','CES'};
SN{end+1} = {'11','YYH'};
SN{end+1} = {'12','LHB'};
SN{end+1} = {'13','JHY'};
SN{end+1} = {'14','KYJ'};
% experiment info. & directory

exp = 'v3';
bhFilePrefix = '_ColorTuningATTEcc_';

nColor = 8; %number of colors presented
nCond = 2; %spatial attention: in vs. out
nTrial = 32; % number of trials within a run
nTarget = 9; % number of target frames within a trial
nCI = 2; % Cardinal vs. Intercardinal
nLMS = 2; % L-M vs. S
carInx = ([1:(nColor/2)]-1)*2+1; %1, 3, 5, 7
intInx = ([1:(nColor/2)])*2; %2,4,6,8
LMInx = [1,5];
SInx = [3,7];

nRun = 8;
nTR = 197;
blank1 = 2;
blank2 = 3;
sTR = 5;
fTR = 1;
tTR = sTR+fTR;
line_colors = {[1 0 0.6], [0.18 0.40 0.73]};

baseDir = sprintf('/Volumes/Duri/data/Color%s', exp);
bhDir = 'BH_data';
csvDir = fullfile(baseDir, 'BH_results');
if ~exist(csvDir, 'dir'); mkdir(csvDir); end

% plot parameters
lineWidth = 3;
lineColors = {[1 0.6 0.784],...
    [1 0 0.6],...
    [0.855 0.702 1],...
    [0.18 0.40 0.73],...
    [0 1 1],...
    [0 0.498 0],...
    [0.4 0.8 0],...
    [0.878 0.537 0.098]};
opacity = 0.3;

Hnames= {'All colors', 'Cardinal colors only'};
CInames = {'Cardinal', 'Int.cardinal'};
LMSnames = {'L-M', 'S'};
attNames = {'In', 'Out'};

attColors = {lineColors{2}, lineColors{4}};
CIColors = {lineColors{2}, lineColors{4}; lineColors{2}, lineColors{4}};
LMSColors = {lineColors{2}, lineColors{4}; lineColors{2}, lineColors{4}};

runNames = {'Run1', 'Run2', 'Run3', 'Run4', 'Run5', 'Run6', 'Run7', 'Run8'};
colorNames = {'Color1', 'Color2', 'Color3', 'Color4', 'Color5', 'Color6', 'Color7', 'Color8'};


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TotalAcc = [];
for xSN = 1:length(SN)
    
    for xRun = 1:nRun
        fileName = strcat(SN{xSN}{1}, SN{xSN}{2}, bhFilePrefix, num2str(xRun), '*.mat');
        fileList = dir(fullfile(baseDir, SN{xSN}{1}, bhDir, fileName));
        if length(fileList) > 1, fprintf('More than two data files for run %d exist!', xRun), end
        load(fullfile(baseDir, SN{xSN}{1}, bhDir, fileList.name));
        targetAppear = saveInfo.RSVP; %frame length. nTrial X nTarget
        
        for xTrial = 1:nTrial
            trialInx = (xTrial-1)*nTarget +1:(xTrial*nTarget);
            trialAcc(xTrial,1) = mean(targetAppear(trialInx, 5));
        end
        condMatrix = saveInfo.Condition; %1column: color, 2column:attention
        condAcc = [ones(nTrial,1)*xSN, ones(nTrial,1)*xRun, condMatrix trialAcc];
        TotalAcc = [TotalAcc; condAcc];
        
        for xColor = 1:nColor
            for xAtt = 1:nCond
                allSubCondAcc{xRun,xColor,xAtt,xSN} =...
                    condAcc(condAcc(:,1) == xColor & condAcc(:,2)==xAtt,3);
            end
        end
        
    end
end

if length(TotalAcc) ~= length(SN)*nRun*nTrial
    error('Some values are missing among runs, trials or subjects!')
end

%%TotalAcc
%%1SN 2Run 3Color 4Attention 5Accuracy

% Run acc.
% for each sub
runAcc = [];
for xSN =1:length(SN)
    for xRun = 1:nRun
        runAcc(xSN,xRun) = mean(TotalAcc(TotalAcc(:,1) == xSN & TotalAcc(:,2) == xRun, 5));
    end
    runSEM(1,1:xRun) = std(runAcc,0,1) ./ sqrt(length(SN)-1);
end

% Color Acc.
colorAcc = [];
for xSN = 1:length(SN)
    for xColor = 1:nColor
        colorAcc(xSN,xColor) = mean(TotalAcc(TotalAcc(:,1) == xSN & TotalAcc(:,3) == xColor, 5));
    end
    colorSEM(1,1:xColor) = std(colorAcc,0,1) ./ sqrt(length(SN)-1);
    
    
end

carColorAcc = mean(colorAcc(:,carInx),2);
carSEM = std(reshape(colorAcc(:,carInx),[],1)) ./ sqrt(length(SN)-1);
intColorAcc = mean(colorAcc(:,intInx),2);
intSEM = std(reshape(colorAcc(:,intInx),[],1)) ./ sqrt(length(SN)-1);
LMColorAcc = mean(colorAcc(:,LMInx),2);
LMSEM = std(reshape(colorAcc(:,LMInx),[],1)) ./ sqrt(length(SN)-1);
SColorAcc = mean(colorAcc(:,SInx),2);
SSEM = std(reshape(colorAcc(:,SInx),[],1)) ./ sqrt(length(SN)-1);

% Attention Acc.
attAcc = [];
for xSN = 1:length(SN)
    for xAtt = 1:nCond
        attAcc(xSN,xAtt) = mean(TotalAcc(TotalAcc(:,1) == xSN & TotalAcc(:,4) == xAtt, 5));
    end
    attSEM(1,1:xAtt) = std(attAcc,0,1) ./ sqrt(length(SN)-1);
end


% Color vs. Attention (2 x 2) acc.
colorAttAcc = [];
for xSN = 1:length(SN)
    for xColor = 1:nColor
        for xAtt = 1:nCond
            colorAttAcc(xSN,xColor,xAtt) =...
                mean(TotalAcc(TotalAcc(:,1) == xSN &...
                TotalAcc(:,3) == xColor &...
                TotalAcc(:,4) == xAtt, 5));
        end
    end
end

CIAttAcc = [];
CIAttAcc(:,1,:) = mean(colorAttAcc(:,carInx,:),2); %bh acc. for Cardinal colors
CIAttAcc(:,2,:) = mean(colorAttAcc(:,intInx,:),2); %bh acc. for Intercardinal colors
LMSAttAcc = [];
LMSAttAcc(:,1,:) = mean(colorAttAcc(:,LMInx,:),2); %bh acc. for L-M colors
LMSAttAcc(:,2,:) = mean(colorAttAcc(:,SInx,:),2); %bh acc. for S colors
attSEM(1,1:xAtt) = std(attAcc,0,1) ./ sqrt(length(SN)-1);


CIAttSEM = []; LMSAttSEM = [];
for xCI = 1:nCI
    for xAtt = 1:nCond
        CIAttSEM(1,xCI,xAtt) = std([CIAttAcc(:,xCI,xAtt)],0,1) / sqrt(length(SN)-1);
        LMSAttSEM(1,xCI,xAtt) = std([LMSAttAcc(:,xCI,xAtt)],0,1) / sqrt(length(SN)-1);
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 1: run
figure(1);
subplot(2, nRun, 1:nRun);
h = plot(mean(runAcc,1), 'color', 'k', 'linewidth', lineWidth);
hold on;
set(gca, 'box', 'off', 'XTickLabel', runNames);
xlim([1 8]); ylim([0 1]);
title(sprintf('%s color blob detection accuracy in each run', exp));
uSEM = []; dSEM = [];
uSEM = mean(runAcc,1)+runSEM(1,:);
dSEM = mean(runAcc,1)-runSEM(1,:);
yAxis = [uSEM fliplr(dSEM)];
xAxis = [1:length(uSEM) fliplr(1:length(dSEM))];
p = patch(xAxis, yAxis, get(h, 'Color'));
set(p, 'FaceAlpha', opacity, 'edgeColor', 'none');
xlim([1 8]); ylim([0 1]);

subplot(2, nColor, nRun+(1:nColor));

% Figure 2: color
b = bar(1:8, mean(colorAcc,1), 'EdgeColor', 'none');
hold on
title(sprintf('%s color blob detection accuracy in each color cond.', exp));
errorbar(1:8, mean(colorAcc,1), colorSEM, colorSEM, 'k', 'linestyle', 'none', 'LineWidth', 1);
set(gca, 'XTickLabel', colorNames, 'Box', 'off');
b.FaceColor = 'flat';
for xColor = 1:nColor
    b.CData(xColor,:) = lineColors{xColor};
end
ylim([0 1]);
hold off
set(gcf, 'renderer', 'painter');


% Figure 2: attention
figure(2);
b = bar(1:2, mean(attAcc,1), 'EdgeColor', 'none');
hold on
title(sprintf('%s color blob detection accuracy\n In vs. Out', exp));
errorbar(1:2, mean(attAcc,1), attSEM, attSEM, 'k', 'linestyle', 'none', 'LineWidth', 1);
set(gca, 'XTickLabel', attNames, 'Box', 'off');
b.FaceColor = 'flat';
for xAtt = 1:nCond
    b.CData(xAtt,:) = attColors{xAtt};
end
ylim([0 1]);
hold off
set(gcf, 'renderer', 'painter');

% Figure 3: color x attention
figure(3);

%Cardinal vs. Intercardinal
subplot(1,2,1);
color(1,:) = mean([CIAttAcc(:,1,1), CIAttAcc(:,1,2)],1); % e.g. cardinal: in,out
color(2,:) = mean([CIAttAcc(:,2,1), CIAttAcc(:,2,2)],1); % e.g. intercardinal: in,out
y = [color(1,:); color(2,:)];
sem(1,:) = [CIAttSEM(1,1,1) CIAttSEM(1,1,2)];
sem(2,:) = [CIAttSEM(1,2,1) CIAttSEM(1,2,2)];
err = [sem(1,:); sem(2,:)];

b = bar(y, 'EdgeColor', 'none');
hold on
for xAtt = 1:nCond
    b(xAtt).FaceColor = CIColors{1,xAtt};
end
set(gca, 'XTickLabel', CInames, 'Box', 'off');
title(sprintf('%s all colors x attention',exp));
[hh, icons, plots, txt] = legend(attNames, 'Location', 'Northeast');
nFactor1 = size(y,1);
nbars = size(y,2);
% Calculate the width for each bar group
barWidth = min(0.8, nbars/(nbars + 1.5));
% Set the position of each error bar in the center of the main bar
for i = 1:nbars
    x = (1:nFactor1) - barWidth/2 + (2*i-1)*barWidth/(2*nbars);
    errorbar(x, y(:,i), err(:,i), err(:,i), 'k', 'linestyle', 'none', 'LineWidth', 1);
end
ylim([0 1]);


%L-M vs. S
subplot(1,2,2);
color(1,:) = mean([LMSAttAcc(:,1,1), LMSAttAcc(:,1,2)],1); % e.g. cardinal: in,out
color(2,:) = mean([LMSAttAcc(:,2,1), LMSAttAcc(:,2,2)],1); % e.g. intercardinal: in,out
y = [color(1,:); color(2,:)];
sem(1,:) = [LMSAttSEM(1,1,1) LMSAttSEM(1,1,2)];
sem(2,:) = [LMSAttSEM(1,2,1) LMSAttSEM(1,2,2)];
err = [sem(1,:); sem(2,:)];

b = bar(y, 'EdgeColor', 'none');
hold on
for xAtt = 1:nCond
    b(xAtt).FaceColor = CIColors{1,xAtt};
end
set(gca, 'XTickLabel', LMSnames, 'Box', 'off');
nFactor1 = size(y,1);
nbars = size(y,2);
% Calculate the width for each bar group
barWidth = min(0.8, nbars/(nbars + 1.5));
% Set the position of each error bar in the center of the main bar
for i = 1:nbars
    x = (1:nFactor1) - barWidth/2 + (2*i-1)*barWidth/(2*nbars);
    errorbar(x, y(:,i), err(:,i), err(:,i), 'k', 'linestyle', 'none', 'LineWidth', 1);
end
ylim([0 1]);
title(sprintf('%s cardinal colors x attention',exp));
set(gcf, 'renderer', 'painter');


% RM ANOVA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('.\n.\n.\n.\n.\n.\n.\n.\n.\n.\n.\n.\n.\n');

% for run (one-way rm anova)
runRMtbl = array2table([[1:length(SN)]', runAcc]);
runRMtbl.Properties.VariableNames = ['SN' runNames];
within = table([1:nRun]', 'VariableNames', {'Run'});
runRMfit = fitrm(runRMtbl, 'Run1-Run8~1', 'WithinDesign', within);
rmANOVA.run = ranova(runRMfit);
if find(0.05 > rmANOVA.run.pValue), sigMsg = 'significant';
elseif find(rmANOVA.run.pValue > 0.05), sigMsg = 'not significant'; end
fprintf('************Run repeated-measures ANOVA results: %s***********\n', sigMsg);
rmANOVA.run

clear sigMsg;

% for colors (one-way rm anova)
colorRMtbl = array2table([[1:length(SN)]', colorAcc]);
colorRMtbl.Properties.VariableNames = ['SN' colorNames];
within = table([1:nColor]', 'VariableNames', {'Color'});
colorRMfit = fitrm(colorRMtbl, 'Color1-Color8~1', 'WithinDesign', within);
rmANOVA.color = ranova(colorRMfit);
if find(0.05 > rmANOVA.color.pValue), sigMsg = 'significant';
elseif find(rmANOVA.color.pValue > 0.05), sigMsg = 'not significant'; end
fprintf('************Color repeated-measures ANOVA results: %s***********\n', sigMsg);
rmANOVA.color
% posthoc.color = multcompare(colorRMfit, 'Color');
% posthoc.color = posthoc.color{find(posthoc.color.Color_1 < posthoc.color.Color_2),:};
fprintf('\n\n')


% for color x attention
% all colors
factorNames = cell(2,2);
factorNames = {'C vs. I', 'Attention'; 'L-M vs. S', 'Attention'};

CIanova(:,1) = reshape(CIAttAcc, [], 1);
CIanova(:,2) = repmat([1:length(SN)]', [nCI*nCond,1]);
CIanova(:,3) = repmat([ones(length(SN),1); ones(length(SN),1)*2], [nCI,1]); % factor1: CI
CIanova(:,4) = [ones(length(SN)*2,1); ones(length(SN)*2,1)*2]; % factor2
rmANOVA.CI = rm_anova2(CIanova(:,1), CIanova(:,2), CIanova(:,3), CIanova(:,4), {factorNames{1,:}});
fprintf('************ %s X attention ANOVA results ************\n', Hnames{1})
rmANOVA.CI
fprintf('\n')

LMSanova(:,1) = reshape(LMSAttAcc, [], 1);
LMSanova(:,2) = repmat([1:length(SN)]', [nLMS*nCond,1]);
LMSanova(:,3) = repmat([ones(length(SN),1); ones(length(SN),1)*2], [nLMS,1]); % factor1: LMS
LMSanova(:,4) = [ones(length(SN)*2,1); ones(length(SN)*2,1)*2]; % factor2
rmANOVA.LMS = rm_anova2(LMSanova(:,1), LMSanova(:,2), LMSanova(:,3), LMSanova(:,4), {factorNames{2,:}});
fprintf('************ %s X attention ANOVA results ************\n', Hnames{2})
rmANOVA.LMS

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save .csv files for JAMOVI...

csvPrefixRun = fullfile(csvDir, sprintf('%s_datRunAcc.csv', exp));
csvPrefixColor = fullfile(csvDir, sprintf('%s_datColorAcc.csv', exp));
csvPrefixCI = fullfile(csvDir, sprintf('%s_datCIAttAcc.csv', exp));
csvPrefixLMS = fullfile(csvDir, sprintf('%s_datLMSAttAcc.csv', exp));

csvRun = fopen(csvPrefixRun, 'a+');
csvColor = fopen(csvPrefixColor, 'a+');
csvCI = fopen(csvPrefixCI, 'a+');
csvLMS = fopen(csvPrefixLMS, 'a+');

% write headers
fprintf(csvRun, 'Run1, Run2, Run3, Run4, Run5, Run6, Run7, Run8\n');
fprintf(csvColor, 'Color1, Color2, Color3, Color4, Color5, Color6, Color7, Color8\n');
fprintf(csvCI, 'C-In, C-Out, I-In, I-Out\n');
fprintf(csvLMS, 'LM-In, S-Out, LM-In, S-Out\n');

% write accuracies
for xSN = 1:length(SN)

fprintf(csvRun, '%2.4f, %2.4f, %2.4f, %2.4f, %2.4f, %2.4f, %2.4f, %2.4f\n',...
    runAcc(xSN,1), runAcc(xSN,2), runAcc(xSN,3), runAcc(xSN,4),...
    runAcc(xSN,5), runAcc(xSN,6), runAcc(xSN,7), runAcc(xSN,8));

fprintf(csvColor, '%2.4f, %2.4f, %2.4f, %2.4f, %2.4f, %2.4f, %2.4f, %2.4f\n',...
    colorAcc(xSN,1), colorAcc(xSN,2), colorAcc(xSN,3), colorAcc(xSN,4),...
    colorAcc(xSN,5), colorAcc(xSN,6), colorAcc(xSN,7), colorAcc(xSN,8));

fprintf(csvCI, '%2.4f, %2.4f, %2.4f, %2.4f\n',...
    CIAttAcc(xSN,1,1), CIAttAcc(xSN,1,2), CIAttAcc(xSN,2,1), CIAttAcc(xSN,2,2)); 

fprintf(csvLMS, '%2.4f, %2.4f, %2.4f, %2.4f\n',...
    LMSAttAcc(xSN,1,1), LMSAttAcc(xSN,1,2), LMSAttAcc(xSN,2,1), LMSAttAcc(xSN,2,2)); 

end




















































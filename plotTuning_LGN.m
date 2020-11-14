% average subjects' data - L-M & S

clear all;
avg_sub={};
avg_sub_card = {};

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

expname = 'v3';
% %
% %RSVP
% SN = {};
% SN{end+1} = '01'; %
% SN{end+1} = '02'; %
% SN{end+1} = '03'; %
% SN{end+1} = '05'; %
% SN{end+1} = '06'; %
% SN{end+1} = '07'; %
% SN{end+1} = '08'; %
% SN{end+1} = '09';
% SN{end+1} = '10'; %
%
% expname = 'RSVP';

nBin = 1;

% root_dir = '/group_hpc/WMShimLab/ColorStudy/v3/';
% root_dir = '/home/sunyoung/psy/Colorv3/';
%root_dir = '/sas2/PECON/PSY/Colorv3/';
% root_dir = '/sas2/PECON/PSY/ColorRSVP/';
root_dir = '/group_hpc/WMShimLab2/PSY_Color/Colorv3/';

addpath('/group_hpc/WMShimLab/PSY_AM_Prediction/analysis_script/');
%root_dir = '/group_hpc/WMShimLab/PSY_ColorStudy/RSVP/';
% root_dir = ['/sas2/PECON/PSY/Color' expname '/'];
fm_dir = 'Img_data/forwardmodel/';
% fm_dir = 'Img_data/forwardmodel_old/';
result_dir='sc_dt_hp_am';
% result_dir='sc_hp';
% result_dir='sc_hp_compareROIs';
% result_name = '_tuning_LGN_.1_roi_JYM_shift.txt';
% result_name = '_tuning_LGN_.1_roi_JYM_scdthp_shift.txt'; % 이게 최종!!!
% result_name = '_tuning_LGN_.1_roi_sc_hp_shift.txt';
% result_name = '_tuning_LGN_.1_shift.txt';
% result_name = 'v3_V1q.001_shift.txt';
% result_name = '_tuning_LGN_.1_JYM_shift.txt';
% result_name = 'v3_V4vq.001_shift.txt';
% result_name = '_tuning_LGN_shift.txt';
% result_name = '_RSVP_V1q.001_shift.txt';
% result_name = '_tuning_LGN_.1_shift.txt';
% result_name = '_RSVP_V4vq.001_shift.txt';
result_name = '_tuning_LGN_hk2_p.05_shift.txt';
% result_name = '_tuning_LGN_hk2_p.05_shift.txt';

f_dir = [fm_dir result_dir];

ROIs = {'LGN'};
nROI = length(ROIs);

nChan = 8;
nCond = 16;

for cc = 1:nCond
    CC{cc} = [];
end

for sub = 1:length(SN)
    
    for cc = 1:nCond
        CC{cc} = [];
    end
    
    clear  TT BB avg_CI card_LM card_S cTT SEM cSEM
    %     fileName = [SN{sub} '_tuning_LGN_shift.txt'];
    fileName = [SN{sub} result_name];
    %     if strcmp(SN{sub}, '04') || strcmp(SN{sub}, '05')
    %     fileName = [SN{sub} '_tuning_LGN_.1_RSVP_shift.txt'];
    %     end
    BB = load(fullfile(root_dir, SN{sub}, f_dir, fileName));
    
    for cc = 1:nCond
        %         CC{cc} = [CC{cc}; (BB(cc,:)+5)./(BB(cc,8)+5) ];
        %         tTT(sub,1:8,cc) = (BB(cc,:)+5)./(BB(cc,8)+5);
        
        CC{cc} = [CC{cc}; zscore(BB(cc,:)) ];
        tTT(sub,1:8,cc) = zscore(BB(cc,:)); %tTT(sub, channel, cond(color*attention))
        
        %         CC{cc} = [CC{cc}; BB(cc,:)-BB(cc,8) ];
        %         tTT(sub,1:8,cc) = BB(cc,:)-BB(cc,8);
    end
    
    % end
    
    % TT(sub,channel,color,cond)
    TT(:,1:8,1,1) = mean(tTT(:,:,[1,3,5,7]),3);
    TT(:,1:8,2,1) = mean(tTT(:,:,[2,4,6,8]),3);
    TT(:,1:8,1,2) = mean(tTT(:,:,[9,11,13,15]),3);
    TT(:,1:8,2,2) = mean(tTT(:,:,[10,12,14,16]),3);
    
    
    for cc = 1:nCond
        Avg(cc,:) = mean(CC{cc},1);
    end
    
    i = 1:4;
    carIdx = (i-1)*2+1;
    intIdx = 2*i;
    avg_CI{1}(1,:) = mean(Avg(carIdx,:),1); %cardinal, in -jyh
    avg_CI{1}(2,:) = mean(Avg(carIdx+8,:),1); %cardinal, out
    avg_CI{2}(1,:) = mean(Avg(intIdx,:),1); %intercardinal, in
    avg_CI{2}(2,:) = mean(Avg(intIdx+8,:),1); %intercardinal, out
    avg_sub_card{sub} = mean([avg_CI{1}(1, :); avg_CI{1}(2,:)]);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % split cardinal into L-M and S
    cTT(:,1:8,1,1) = mean(tTT(:,:,[1,5]),3);
    cTT(:,1:8,2,1) = mean(tTT(:,:,[3,7]),3);
    cTT(:,1:8,1,2) = mean(tTT(:,:,[9,13]),3);
    cTT(:,1:8,2,2) = mean(tTT(:,:,[11,15]),3);
    
    LMIdx = [1 5];
    SIdx = [3 7];
    card_LM(1,:) = mean(Avg(LMIdx,:),1);
    card_LM(2,:) = mean(Avg(LMIdx+8,:),1);
    card_S(1,:) = mean(Avg(SIdx,:),1);
    card_S(2,:) = mean(Avg(SIdx+8,:),1);
    
    cSEM = std(cTT,0,1)./sqrt(length(SN)-1);
    
    % add ch 8 before ch 8
    card_LM(:,2:9) = card_LM(:,1:8);
    card_LM(:,1) = card_LM(:,9);
    card_S(:,2:9) = card_S(:,1:8);
    card_S(:,1) = card_S(:,9);
    card{1} = card_LM;
    card{2} = card_S;
    for col = 1:2
        cSEM(1,2:9,col,1) =  cSEM(1,1:8,col,1);
        cSEM(1,1,col,1) =  cSEM(1,9,col,1);
        cSEM(1,2:9,col,2) =  cSEM(1,1:8,col,2);
        cSEM(1,1,col,2) =  cSEM(1,9,col,2);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % % plot % %
    
    % when shifted, tuning functions are centered on channel 4
    % linecolors = {'r', '--r','b', '--b', 'c', '--c','m', '--m'};
    linecolors = {[1 0.6 0.784], [1 0 0.6], [0.855 0.702 1], [0.18 0.40 0.73], ...
        [0 1 1], [0 0.498 0], [0.4 0.8 0], [0.878 0.537 0.098]};
    
    %calculate SEM
    SEM = std(TT,0,1)./sqrt(length(SN)-1);
    
    for col = 1:2
        % add ch 8 before ch 8
        avg_CI{col}(:,2:9) = avg_CI{col}(:,1:8);
        avg_CI{col}(:,1) = avg_CI{col}(:,9);
        SEM(1,2:9,col,1) =  SEM(1,1:8,col,1);
        SEM(1,1,col,1) =  SEM(1,9,col,1);
        SEM(1,2:9,col,2) =  SEM(1,1:8,col,2);
        SEM(1,1,col,2) =  SEM(1,9,col,2);
    end
    avg_sub{sub} = mean([avg_CI{1}; avg_CI{2}]);
    
    %sub_graph(xsn,:) = mean([avg_CI{1}; avg_CI{2}]);
    
    
    figure(1); %clf;
    avg_colors = {linecolors{2}, linecolors{4}};
    color_names = {'Cardinal', 'Intercardinal'};
    opacity = 0.3;
    expcond = 2;
    colcond = 2;
    for col = 1:colcond
        % subplot(2,2,col);
        
        subplot(4,length(SN), (col-1)*length(SN)+sub)
        for cc = 1:expcond
            
            TempUpperPart=[]; UpperPart=[]; TempLowerPart=[]; LowerPart=[]; XAxis=[]; YAxis=[];
            TempUpperPart=avg_CI{col}(cc,:)+SEM(1,:,col,cc);
            UpperPart=TempUpperPart(1:length(find(~isnan(TempUpperPart))));
            TempLowerPart=avg_CI{col}(cc,:)-SEM(1,:,col,cc);
            LowerPart=TempLowerPart(1:length(find(~isnan(TempLowerPart))));
            YAxis=[UpperPart fliplr(LowerPart)];
            TimeRange=1:length(UpperPart);
            XAxis=[TimeRange fliplr(TimeRange)];
            PlotHandle=plot(avg_CI{col}(cc,:), 'color', avg_colors{cc},'linewidth', 3); hold on
            % PatchHandle=patch(XAxis,YAxis,get(PlotHandle,'Color'));
            % set(PatchHandle,'FaceAlpha',opacity);
            % set(PatchHandle,'EdgeColor','none');
        end
        
        xlim([0 10]);
        ylim([-2 2]);
        % ylim([0.96 1.08]);
        if col == 1
            title([ 's' SN{sub} ' ' color_names{col}]);
        elseif col == 2
            title(color_names{col});
        end
        % legend({'Inner', '', 'Outer'}, 'location', 'NorthEast');
        
    end
    % suptitle(['subj' SN{sub}]);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % figure(sub); %clf;
    avg_colors = {linecolors{2}, linecolors{4}};
    color_names = {'L-M', 'S'};
    opacity = 0.3;
    expcond = 2;
    colcond = 2;
    for col = 1:colcond
        % subplot(2,2,col+2);
        
        subplot(4,length(SN), (col+1)*length(SN)+sub)
        for cc = 1:expcond
            TempUpperPart=[]; UpperPart=[]; TempLowerPart=[]; LowerPart=[]; XAxis=[]; YAxis=[];
            TempUpperPart=card{col}(cc,:)+cSEM(1,:,col,cc);
            UpperPart=TempUpperPart(1:length(find(~isnan(TempUpperPart))));
            TempLowerPart=card{col}(cc,:)-cSEM(1,:,col,cc);
            LowerPart=TempLowerPart(1:length(find(~isnan(TempLowerPart))));
            YAxis=[UpperPart fliplr(LowerPart)];
            TimeRange=1:length(UpperPart);
            XAxis=[TimeRange fliplr(TimeRange)];
            PlotHandle=plot(card{col}(cc,:), 'color', avg_colors{cc},'linewidth', 3); hold on
            % PatchHandle=patch(XAxis,YAxis,get(PlotHandle,'Color'));
            % set(PatchHandle,'FaceAlpha',opacity);
            % set(PatchHandle,'EdgeColor','none');
        end
        
        xlim([0 10]);
        ylim([-2 2]);
        % ylim([0.96 1.08]);
        title(color_names{col});
        % legend({'Inner', '', 'Outer'}, 'location', 'NorthEast');
        hold off;
        
    end
end

set(gcf,'renderer','Painters');
%%

% average subjects' data - L-M & S

% root_dir = '/group_hpc/WMShimLab/ColorStudy/v3/';
% root_dir = '/home/sunyoung/psy/Colorv3/';
% root_dir = '/sas2/PECON/PSY/Colorv3/';
% root_dir = '/sas2/PECON/PSY/ColorRSVP/';
% fm_dir = 'Img_data/forwardmodel/';
% result_dir='sc_dt_hp_am';
f_dir = [fm_dir result_dir];

ROIs = {'LGN'};
nROI = length(ROIs);

nChan = 8;
nCond = 16;

for cc = 1:nCond
    CC{cc} = [];
end
clear tTT TT BB avg_CI card_LM card_S cTT SEM cSEM
for sub = 1:length(SN)
    %     fileName = [SN{sub} '_tuning_LGN_shift.txt'];
    fileName = [SN{sub} result_name];
    %     if strcmp(SN{sub}, '04') || strcmp(SN{sub}, '05')
    %     fileName = [SN{sub} '_tuning_LGN_.1_RSVP_shift.txt'];
    %     end
    %     fileName = [SN{sub} '_tuning_LGN_.1_RSVP_shift.txt'];
    BB = load(fullfile(root_dir, SN{sub}, f_dir, fileName));
    
    for cc = 1:nCond
        %         CC{cc} = [CC{cc}; (BB(cc,:)+5)./(BB(cc,8)+5) ];
        %         tTT(sub,1:8,cc) = (BB(cc,:)+5)./(BB(cc,8)+5);
        
        CC{cc} = [CC{cc}; zscore(BB(cc,:)) ];
        tTT(sub,1:8,cc) = zscore(BB(cc,:));
        %
        %         CC{cc} = [CC{cc}; BB(cc,:)-BB(cc,8) ];
        %         tTT(sub,1:8,cc) = BB(cc,:)-BB(cc,8);
    end
end

% TT(sub,channel,color,cond)
TT(:,1:8,1,1) = mean(tTT(:,:,[1,3,5,7]),3);
TT(:,1:8,2,1) = mean(tTT(:,:,[2,4,6,8]),3);
TT(:,1:8,1,2) = mean(tTT(:,:,[9,11,13,15]),3);
TT(:,1:8,2,2) = mean(tTT(:,:,[10,12,14,16]),3);

for cc = 1:nCond
    Avg(cc,:) = mean(CC{cc},1);
end

i = 1:4;
carIdx = (i-1)*2+1;
intIdx = 2*i;
avg_CI{1}(1,:) = mean(Avg(carIdx,:),1);
avg_CI{1}(2,:) = mean(Avg(carIdx+8,:),1);
avg_CI{2}(1,:) = mean(Avg(intIdx,:),1);
avg_CI{2}(2,:) = mean(Avg(intIdx+8,:),1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% split cardinal into L-M and S
cTT(:,1:8,1,1) = mean(tTT(:,:,[1,5]),3);
cTT(:,1:8,2,1) = mean(tTT(:,:,[3,7]),3);
cTT(:,1:8,1,2) = mean(tTT(:,:,[9,13]),3);
cTT(:,1:8,2,2) = mean(tTT(:,:,[11,15]),3);

LMIdx = [1 5];
SIdx = [3 7];
card_LM(1,:) = mean(Avg(LMIdx,:),1);
card_LM(2,:) = mean(Avg(LMIdx+8,:),1);
card_S(1,:) = mean(Avg(SIdx,:),1);
card_S(2,:) = mean(Avg(SIdx+8,:),1);

cSEM = std(cTT,0,1)./sqrt(length(SN)-1);

% add ch 8 before ch 8
card_LM(:,2:9) = card_LM(:,1:8);
card_LM(:,1) = card_LM(:,9);
card_S(:,2:9) = card_S(:,1:8);
card_S(:,1) = card_S(:,9);
card{1} = card_LM;
card{2} = card_S;
for col = 1:2
    cSEM(1,2:9,col,1) =  cSEM(1,1:8,col,1);
    cSEM(1,1,col,1) =  cSEM(1,9,col,1);
    cSEM(1,2:9,col,2) =  cSEM(1,1:8,col,2);
    cSEM(1,1,col,2) =  cSEM(1,9,col,2);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % plot % %

% when shifted, tuning functions are centered on channel 4
% linecolors = {'r', '--r','b', '--b', 'c', '--c','m', '--m'};
linecolors = {[1 0.6 0.784], [1 0 0.6], [0.855 0.702 1], [0.18 0.40 0.73], ...
    [0 1 1], [0 0.498 0], [0.4 0.8 0], [0.878 0.537 0.098]};

%calculate SEM
SEM = std(TT,0,1)./sqrt(length(SN)-1);

for col = 1:2
    % add ch 8 before ch 8
    avg_CI{col}(:,2:9) = avg_CI{col}(:,1:8);
    avg_CI{col}(:,1) = avg_CI{col}(:,9);
    SEM(1,2:9,col,1) =  SEM(1,1:8,col,1);
    SEM(1,1,col,1) =  SEM(1,9,col,1);
    SEM(1,2:9,col,2) =  SEM(1,1:8,col,2);
    SEM(1,1,col,2) =  SEM(1,9,col,2);
end

figure(); clf;
avg_colors = {linecolors{2}, linecolors{4}};
color_names = {'Cardinal', 'Intercardinal'};
opacity = 0.3;
expcond = 2;
colcond = 2;
for col = 1:colcond
    subplot(2,2,col);
    
    % subplot(4,length(SN), (col-1)*length(SN)+sub)
    for cc = 1:expcond
        
        TempUpperPart=[]; UpperPart=[]; TempLowerPart=[]; LowerPart=[]; XAxis=[]; YAxis=[];
        TempUpperPart=avg_CI{col}(cc,:)+SEM(1,:,col,cc);
        UpperPart=TempUpperPart(1:length(find(~isnan(TempUpperPart))));
        TempLowerPart=avg_CI{col}(cc,:)-SEM(1,:,col,cc);
        LowerPart=TempLowerPart(1:length(find(~isnan(TempLowerPart))));
        YAxis=[UpperPart fliplr(LowerPart)];
        TimeRange=1:length(UpperPart);
        XAxis=[TimeRange fliplr(TimeRange)];
        PlotHandle=plot(avg_CI{col}(cc,:), 'color', avg_colors{cc},'linewidth', 3); hold on
        PatchHandle=patch(XAxis,YAxis,get(PlotHandle,'Color'));
        set(PatchHandle,'FaceAlpha',opacity);
        set(PatchHandle,'EdgeColor','none');
    end
    
    xlim([0 10]);
    ylim([-1 1]);
    % ylim([-2 2]);
    % ylim([0.96 1.08]);
    % if col == 1
    % title([ 's' SN{sub} ' ' color_names{col}]);
    % elseif col == 2
    title(color_names{col});
    % end
    % legend({'Inner', '', 'Outer'}, 'location', 'NorthEast');
    
end
% suptitle(['subj' SN{sub}]);

%cSEM_sub = std(subGraph(xsn,:),0,1)./sqrt(length(SN)-1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure(sub); %clf;
avg_colors = {linecolors{2}, linecolors{4}};
color_names = {'L-M', 'S'};
opacity = 0.3;
expcond = 2;
colcond = 2;
for col = 1:colcond
    subplot(2,2,col+2);
    
    % subplot(4,length(SN), (col+1)*length(SN)+sub)
    for cc = 1:expcond
        TempUpperPart=[]; UpperPart=[]; TempLowerPart=[]; LowerPart=[]; XAxis=[]; YAxis=[];
        TempUpperPart=card{col}(cc,:)+cSEM(1,:,col,cc);
        UpperPart=TempUpperPart(1:length(find(~isnan(TempUpperPart))));
        TempLowerPart=card{col}(cc,:)-cSEM(1,:,col,cc);
        LowerPart=TempLowerPart(1:length(find(~isnan(TempLowerPart))));
        YAxis=[UpperPart fliplr(LowerPart)];
        TimeRange=1:length(UpperPart);
        XAxis=[TimeRange fliplr(TimeRange)];
        PlotHandle=plot(card{col}(cc,:), 'color', avg_colors{cc},'linewidth', 3); hold on
        PatchHandle=patch(XAxis,YAxis,get(PlotHandle,'Color'));
        set(PatchHandle,'FaceAlpha',opacity);
        set(PatchHandle,'EdgeColor','none');
    end
    
    xlim([0 10]);
    ylim([-1 1]);
    % ylim([0.96 1.08]);
    title(color_names{col});
    % legend({'Inner', '', 'Outer'}, 'location', 'NorthEast');
    hold off;
    
end

% % quantify by correlation with channel basis function
% x = linspace(0,pi,9);
% pred = sin(x).^6;
% pred = pred./max(pred);
% % plot(pred)
%
% % TT(sub,channel,color,cond)
% xCorr = []; xP=[];
% for color = 1:2
%     for cond = 1:2
%         for sub = 1:length(SN)
%         [xCorr(sub,color,cond) xP(sub,color,cond)] = corr(pred',[TT(sub,8,color,cond) TT(sub,:,color,cond)]');
%
%         end
%     end
% end

% fprintf('===Correlation with centered basis function===\n');
% fprintf('car-in\tcar-out\tint-in\tint-out\n');
%
% for color =1:2
%     for cond = 1:2
%         fprintf('%1.3f\t', mean(xCorr(:,color,cond),1));
%     end
% end
% fprintf('\n\n');

% information fidelity %
%     % FIDELITY %
% clear xSlope data e
xFid = [];
for color = 1:2
    for cond = 1:2
        data = TT(:,:,color,cond);
        %
        nChan = size(data,2);
        nSub = size(data,1);
        c_center = ceil(nChan/2); % for shifted data
        %
        ori_unit = 2*pi/nChan;
        ori_rad = 0:ori_unit:2*pi;
        ori_rad = ori_rad-ori_rad(c_center);
        e = zeros(nSub,nChan);
        for i = 1:nChan
            e(:,i)=data(:,i).*cos(abs(ori_rad(i)));
        end
        %
        xFid(:,color,cond) = mean(e,2);
    end
end
% shape data for anova
colcond=2; nCond = 2;
xAnova = reshape(xFid,[length(SN)*colcond*nCond 1]);
xAnova(:,2) = repmat((1:length(SN))',[colcond*nCond 1]);
temp = mod(((1:length(xAnova))-1),length(xAnova)/2)+1;
xAnova(:,3) = ceil(temp./(length(xAnova)/nCond/colcond));
xAnova(:,4) = ceil((1:length(xAnova))./(length(xAnova)/2));
factornames = {'Color', 'Attention'};

fprintf('===ANOVA all color fidelity===\n');
AnovaTable = rm_anova2(xAnova(:,1),xAnova(:,2),xAnova(:,3),xAnova(:,4),factornames)
An_cycle_f(i) = AnovaTable{2,5};
An_cycle_p(i) = AnovaTable{2,6};
An_cond_f(i) = AnovaTable{3,5};
An_cond_p(i) = AnovaTable{3,6};
An_inter_f(i) = AnovaTable{4,5};
An_inter_p(i) = AnovaTable{4,6};


% information fidelity %
%        FIDELITY %
clear xSlope data e
xFid = [];
for color = 1:2
    for cond = 1:2
        data = cTT(:,:,color,cond);
        
        nChan = size(data,2);
        nSub = size(data,1);
        c_center = ceil(nChan/2); % for shifted data
        
        ori_unit = 2*pi/nChan;
        ori_rad = 0:ori_unit:2*pi;
        ori_rad = ori_rad-ori_rad(c_center);
        e = zeros(nSub,nChan);
        for i = 1:nChan
            e(:,i)=data(:,i).*cos(abs(ori_rad(i)));
        end
        
        xFid(:,color,cond) = mean(e,2);
    end
end
% shape data for anova
colcond=2; nCond = 2;
xAnova = reshape(xFid,[length(SN)*colcond*nCond 1]);
xAnova(:,2) = repmat((1:length(SN))',[colcond*nCond 1]);
temp = mod(((1:length(xAnova))-1),length(xAnova)/2)+1;
xAnova(:,3) = ceil(temp./(length(xAnova)/nCond/colcond));
xAnova(:,4) = ceil((1:length(xAnova))./(length(xAnova)/2));
factornames = {'Color', 'Attention'};
fprintf('===ANOVA car color fidelity===\n');
AnovaTable = rm_anova2(xAnova(:,1),xAnova(:,2),xAnova(:,3),xAnova(:,4),factornames)
An_cycle_f(i) = AnovaTable{2,5};
An_cycle_p(i) = AnovaTable{2,6};
An_cond_f(i) = AnovaTable{3,5};
An_cond_p(i) = AnovaTable{3,6};
An_inter_f(i) = AnovaTable{4,5};
An_inter_p(i) = AnovaTable{4,6};

% linear regression %
% TT(sub,channel,color,cond)
% centered on 4
folded = zeros(size(TT,1), 5, size(TT,3), size(TT,4));
folded(:,1,:,:) = TT(:,8,:,:);
folded(:,2,:,:) = (TT(:,1,:,:)+TT(:,7,:,:))/2;
folded(:,3,:,:) = (TT(:,2,:,:)+TT(:,6,:,:))/2;
folded(:,4,:,:) = (TT(:,3,:,:)+TT(:,5,:,:))/2;
folded(:,5,:,:) = TT(:,4,:,:);
x = 1:size(folded,2);
LinFit = zeros(size(TT,1), 2, size(TT,3), size(TT,4));
for color = 1:2
    for cond = 1:2
        for sub = 1:size(folded,1)
            y=folded(sub,:,color,cond);
            P = polyfit(x,y,1);
            %             rsq(sub,(color-1)*2+cond) = linrsq(P,x,y);
            LinFit(sub,:,color,cond) = P;
        end
        %         CC=mean(LinFit)
    end
end

% bucket slope
for bb = 1:2
    for cc = 1:2
        xSlope(:,cc,bb) = LinFit(:,1,cc,bb);
    end
end

xSlope_all = xSlope;

% shape data for anova
colcond=2; nCond = 2;
xAnova = reshape(xSlope,[length(SN)*colcond*nCond 1]);
xAnova(:,2) = repmat((1:length(SN))',[colcond*nCond 1]);
temp = mod(((1:length(xAnova))-1),length(xAnova)/2)+1;
xAnova(:,3) = ceil(temp./(length(xAnova)/nCond/colcond));
xAnova(:,4) = ceil((1:length(xAnova))./(length(xAnova)/2));
factornames = {'Color', 'Attention'};
fprintf('===ANOVA all color slope===\n');
AnovaTable = rm_anova2(xAnova(:,1),xAnova(:,2),xAnova(:,3),xAnova(:,4),factornames)


% %t-test with 0
% for color = 1:2
%     for cond = 1:2
% [H, p((color-1)*2+cond)] = ttest(LinFit(:,1,color,cond), ...
%     zeros(length(LinFit(:,1,color,cond)),1), 'tail','both');
%     end
% end
%
% % t-test attention effect
% [H, P(1)] = ttest(LinFit(:,1,1,1), LinFit(:,1,1,2), 'tail','both');
% [H, P(2)] = ttest(LinFit(:,1,2,1), LinFit(:,1,2,2), 'tail','both');
%
% fprintf('Stats n = %d',length(SN));
% fprintf('\n==Significant tuning: %s==\n',result_dir);
% fprintf('car-in\tcar-out\tint-in\tint-out\n');
% fprintf('%1.3f\t%1.3f\t%1.3f\t%1.3f\n\n', p(1), p(2), p(3), p(4));
% fprintf('====Attention Effect====\n');
% fprintf('Cardinal\tIntercardinal\n');
% fprintf('%1.3f\t\t%1.3f\n\n', P(1), P(2));

TT=cTT;
folded = zeros(size(TT,1), 5, size(TT,3), size(TT,4));
folded(:,1,:,:) = TT(:,8,:,:);
folded(:,2,:,:) = (TT(:,1,:,:)+TT(:,7,:,:))/2;
folded(:,3,:,:) = (TT(:,2,:,:)+TT(:,6,:,:))/2;
folded(:,4,:,:) = (TT(:,3,:,:)+TT(:,5,:,:))/2;
folded(:,5,:,:) = TT(:,4,:,:);
x = 1:size(folded,2);
LinFit = zeros(size(TT,1), 2, size(TT,3), size(TT,4));
for color = 1:2
    for cond = 1:2
        for sub = 1:size(folded,1)
            y=folded(sub,:,color,cond);
            P = polyfit(x,y,1);
            %             rsq(sub,(color-1)*2+cond) = linrsq(P,x,y);
            LinFit(sub,:,color,cond) = P;
        end
        %         CC=mean(LinFit)
    end
end


% bucket slope
for bb = 1:2
    for cc = 1:2
        xSlope(:,cc,bb) = LinFit(:,1,cc,bb);
    end
end

xSlope_car = xSlope;

% shape data for anova
colcond=2; nCond = 2;
xAnova = reshape(xSlope,[length(SN)*colcond*nCond 1]);
xAnova(:,2) = repmat((1:length(SN))',[colcond*nCond 1]);
temp = mod(((1:length(xAnova))-1),length(xAnova)/2)+1;
xAnova(:,3) = ceil(temp./(length(xAnova)/nCond/colcond));
xAnova(:,4) = ceil((1:length(xAnova))./(length(xAnova)/2));
factornames = {'Color', 'Attention'};
fprintf('===ANOVA car color slope===\n');
AnovaTable = rm_anova2(xAnova(:,1),xAnova(:,2),xAnova(:,3),xAnova(:,4),factornames)
%
% %t-test with 0
% for color = 1:2
%     for cond = 1:2
% [H, p((color-1)*2+cond)] = ttest(LinFit(:,1,color,cond), ...
%     zeros(length(LinFit(:,1,color,cond)),1), 'tail','both');
%     end
% end
%
% % t-test attention effect
% [H, P(1)] = ttest(LinFit(:,1,1,1), LinFit(:,1,1,2), 'tail','both');
% [H, P(2)] = ttest(LinFit(:,1,2,1), LinFit(:,1,2,2), 'tail','both');
%
% fprintf('\n==Significant tuning: %s==\n',result_dir);
% fprintf('LM-in\tLM-out\tS-in\tS-out\n');
% fprintf('%1.3f\t%1.3f\t%1.3f\t%1.3f\n\n', p(1), p(2), p(3), p(4));
% fprintf('====Attention Effect====\n');
% fprintf('LM\t\tS\n');
% fprintf('%1.3f\t\t%1.3f\n\n', P(1), P(2));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






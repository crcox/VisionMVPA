% basic setup
% cd /Volumes/LKCLAB/MRI/Motion3D/Adaptation
% cd /Volumes/LKCLAB/MRI/Raw
cd /home/chris/VisionMVPA
addpath('/home/chris/VisionMVPA');

% cd Y:\MRI\Raw
close all
% clear all

% analysis params
normFlag = 0;
cothresh = 0.3; % threshold for average voxel activity across all conditions % could be 0 for example
phWindow = [0 2*pi];%pi/2]; %phase window - can go from 0-2pi % gets used in line 92 when doing the inplane analysis

% subtype = 'tmscr';
% EXTENDED code so it runs all conditions at the same time
% s = {'053008c','053008c'}; % cormack, twice for bug
% cormack, huk, rokers
% s = {'060208a','060208b','060208c',... % DRDS wedges and anti-correlated
%     '061008a','061008b','061008c',... % DRDS wedges and spatial scramble
%     '061308a','061308b','061308c'}; % DRDS wedges and temporal scramble
s = {'Pilot-BR-002-Take-2'}; %{'121108b'}; %s = {'121108a'};


% sub-index in to do single-subject analyses
si = 1:length(s); %for all subjects
s = s(si);

nConds = 3;

c{1} = {[1:10]}; %{[1 2]}; % which scans to use
nScans = 10;%10; 

% rois = {'V1','MT+','IPSm'};
% rois = {'V1','V2','V3','V3A','V4v','LO','MT+','IPSp','IPSm','IPSa'};
% rois = {'V1','V2','V3','V3A','V4v','LO','MT+','MT','MST'};
rois = {'V1'}; %{'MT+'};


fName = ['perceptualadaptation_' s{1} '_' num2str(cothresh*100) '_' rois{1}];
datadir = '/home/chris/VisionMVPA';

ok = 1; % error-checking flag
for k=1 %:4 % Number of conditions
    ts{k} = deal(cell(length(rois),1));
end


%% loop thru sessions and do analysis
for i = 1:length(s)
    
    % go to each session and init data structures
    disp(s{i});
    chdir(fullfile(datadir,s{i}));
    inplane = [];
    
    for j = 1:length(rois)
        
        % Compute average phase of MT localizer
        inplane = initHiddenInplane('Averages'); % Done so that you only select voxels that are active across all scans on average
        refScan = 1; % Refscan is average of scan 1 and 10
        
        % load ROI
        [inplane, ok] = loadROI(inplane, rois{j});
        if ~ok
            error('ROI loading failed!')
        end
        
        inplane = loadCorAnal(inplane);
        inplane = restrictROI(inplane,refScan,cothresh,phWindow);
        
        numberOfUseableVoxels = length(inplane.ROIs.coords),
        
        % Switch to motioncomp scans
        inplane = selectDataType(inplane,'MotionComp');
        inplane = loadCorAnal(inplane);
        
        for k = 1:nScans % Number of scans/session
            tSeries{k} = meanTSeries(inplane,k);
        end
        
        % now stuff these into appropriate data structures
        disp(['Analyzing ' rois{j} ' data']);
        
        for k = 1 %:4 % number of conditions
            ts{k}{j} = [ts{k}{j} tSeries(c{k}{i})];
        end
    end
end
disp('Extracted raw data')

%% Reshape time series
% These should be converted to some kind of viewGet() calls.
nCycles   = 36; % 36 - trials in experiment;
nRepeats = 1; % # repeats of same stimulus / session
nFrames   = 180; % 175;
framesPerCycle = 15; % number of frames we want to see in the reshaped time series 180/36 = %15; %nFrames/nCycles;


nScans = length(ts{k}{i});


ts = EventRelated(rois,nScans,ts);

% our experiment
d = dir('Behavior/*.mat');
% for i = 1:length(rois)
%     for j = 1:nScans % 10
%         load(['Behavior/' d(c{1}{1}(j)).name]); % c{1}{1}-1 % NEEDS TO BE FIXED (parametrized)
%         congruentTs{1}{i}{j} = ts{1}{i}{j}(find(pa.response(3:end-2,3)==1),:); % 0),:); %opposite probe
%         congruentTs{2}{i}{j} = ts{1}{i}{j}(find(pa.response(3:end-2,3)==2),:); %-1),:); %same probe
%         if nConds > 2
%             congruentTs{3}{i}{j} = ts{1}{i}{j}(find(pa.response(3:end-2,3)==3),:); % blank probe
%         end
%     end
% end
for i = 1:length(rois)
    for j = 1:nScans % 10
        load(['Behavior/' d(c{1}{1}(j)).name]); % c{1}{1}-1 % NEEDS TO BE FIXED (parametrized)
        congruentTs{1}{i}{j} = ts{1}{i}{j}(find(pa.response(1:end-2,3)==1),:); % 0),:); %opposite probe
        congruentTs{2}{i}{j} = ts{1}{i}{j}(find(pa.response(1:end-2,3)==2),:); %-1),:); %same probe
        if nConds > 2
            congruentTs{3}{i}{j} = ts{1}{i}{j}(find(pa.response(1:end-2,3)==3),:); % blank probe
        end
    end
end


for k = 1:nConds %:4 % # conditions
    for i = 1:length(rois)
        % Get number of scans in this condition
        for j = 1:nScans % Number of sessions / condition
            %Compute the average single cycle
            singleCycleTs{k}{i}{j} = mean(congruentTs{k}{i}{j});
            
            howMany = size(congruentTs{k}{i}{j});
            singleCycleStdErrTs{k}{i}{j} = std(congruentTs{k}{i}{j})/(sqrt(howMany(1)));
        end
        % framesPerCycle=framesPerCycle+1;
        rsSingleCycleTs{k}{i} = reshape(cell2mat(singleCycleTs{k}{i}),[framesPerCycle nScans]);
        avgSingleCycleTs{k}{i} = mean(rsSingleCycleTs{k}{i}'); % One of these is framesPerCycle+1, the other number of condition repeats
        avgSingleCycleStdErrTs{k}{i} = std(rsSingleCycleTs{k}{i}')./sqrt(nScans*length(congruentTs{k}{i})); %(nScans); 
        % THIS WAS A VERY CONSERVATIVE ESTIMATE OF STD ERR (based on nScans,
        % rather than total trials (nScans x nTrials)
    end
end

% aveAllCondsTs = mean([rsSingleCycleTs{1}{1}; rsSingleCycleTs{2}{1}; rsSingleCycleTs{3}{1}]); %mean([avgSingleCycleTs{1}{1}; avgSingleCycleTs{2}{1}; avgSingleCycleTs{3}{1}]);
% stderrAllCondsTs = mean([rsSingleCycleTs{1}{1}; rsSingleCycleTs{2}{1}; rsSingleCycleTs{3}{1}])/sqrt(nScans*(length(rsSingleCycleTs{1}{1})+length(rsSingleCycleTs{2}{1})+length(rsSingleCycleTs{3}{1})));  %mean([avgSingleCycleTs{1}{1}; avgSingleCycleTs{2}{1}; avgSingleCycleTs{3}{1}])/sqrt(3);

% aveAllStimCondsTs = mean([rsSingleCycleTs{1}{1}'; rsSingleCycleTs{2}{1}']); %mean([avgSingleCycleTs{1}{1}; avgSingleCycleTs{2}{1}; avgSingleCycleTs{3}{1}]);
% stderrAllStimCondsTs = mean([rsSingleCycleTs{1}{1}'; rsSingleCycleTs{2}{1}'])/sqrt(nScans*(length(rsSingleCycleTs{1}{1})+length(rsSingleCycleTs{2}{1})));  %mean([avgSingleCycleTs{1}{1}; avgSingleCycleTs{2}{1}; avgSingleCycleTs{3}{1}])/sqrt(3);



disp('Built timeseries')

%% save out results
outDir = '/home/chris/VisionMVPA/Outputs';
if ~exist(outDir,7)
    mkdirp(outDir);
end
save(fullfile(outDir,fName)); % save the whole damn workspace!

%% load results
% if ~exist('cotresh','var')
%     cothresh = 0.3;
% end
% WARNING: This will change the file you intend to load
% fName = ['perceptualadaptation_' s{1} '_' num2str(cothresh*100) '_' rois{1}];
load(fullfile(outDir,fName));

%% plot results (Baseline subtraction)

figDir = '/home/chris/VisionMVPA/Figures';
if ~exist(figDir,7)
    mkdirp(figDir);
end
close all;

if nConds > 2
    
    LineWidth = 1.5;
    MarkerSize = 8; %10;
    FontSize = 14;
    timeList = 0:1.5:1.5*14; %-10*1.5:1.5:1.5*14; %0:1.5:1.5*14; %14;
    series_color = {[0 1 0], [1 0 0],[0 0 1]};
    
    series_style = {'-','--','---'};
    
    for j = 1:length(rois)
        %     subplot(3,3,j) % no subplots for one ROI
        for k = 1:2 %these are conditions
            hold on
            eh = errorbar(timeList,avgSingleCycleTs{k}{j}-avgSingleCycleTs{3}{j},avgSingleCycleStdErrTs{k}{j},series_style{k},'Color',[0 0 0],'LineWidth',LineWidth);
            plot([timeList(1) timeList(end)], [0 0], 'k--')
            %         removeErrorBarEnds(eh);
            h{k} = plot(timeList, avgSingleCycleTs{k}{j}-avgSingleCycleTs{3}{j},'o','MarkerEdgeColor',0*[1 1 1]);
            
            set(h{k},'Marker','o','MarkerSize',MarkerSize,'MarkerFaceColor',series_color{k},'LineWidth',LineWidth);
            set(gca,'LineWidth',LineWidth)

            set(gca,'YTick',[-2:.1:2]);
            set(gca,'XTick',[0 7.5 15 22.5]);
            set(gca,'XLim',[-2 17.5]);
            set(gca,'FontSize',FontSize)
            
        end
        legend([h{1}(1) h{2}(1)],'opposite','same','Location','SouthWest')
        legend boxoff;
    end
    box on;
    xlabel('Time (s)','FontSize',24)
    ylabel('Relative fMRI response (% change in BOLD)','FontSize',24)
    set(gcf,'Color','w');
    set(gcf,'Position',[100 100 750 500]);
    orient landscape
    set(gcf,'PaperPosition',[0 0 11 8.5]);
    th = title([rois{1} 'Adaptation subtracted. Threshold ='  num2str(cothresh) 'Phase window = ' num2str(phWindow(2)*(180/pi))]);
    set(th,'FontAngle','italic');
    set(th,'FontWeight','bold');
%     saveas(gcf,fullfile(figDir,['perceptualadaptation_subtract',s{1},'_',num2str(cothresh*100),'_',rois{1}]),'fig');
    saveas(gcf,fullfile(figDir,['perceptualadaptation_subtract',s{1},'_',num2str(cothresh*100),'_',rois{1}]),'pdf');
end

%% plot results (No baseline subtraction)
% close all;
LineWidth = 1.5; %1.5;
MarkerSize = 8; %8; %10;
FontSize = 14; %14;
series_color = {[0 1 0], [1 0 0],[0 0 1],[1 1 0]};

series_style = {'-','--','--','-'};
figure, 
hold on
yh = ylabel('fMRI response (% BOLD)');
for j = 1:length(rois)
    %figure, hold on
    %     subplot(3,3,j)
    for k = 1:nConds %[1 3] %these are conditions
        hold on
        eh = errorbar(timeList,avgSingleCycleTs{k}{j},avgSingleCycleStdErrTs{k}{j},series_style{k},'Color',[0 0 0],'LineWidth',LineWidth);
        plot([timeList(1) timeList(end)], [0 0], 'k--')

        h{k} = plot(timeList, avgSingleCycleTs{k}{j},'o','MarkerEdgeColor',0*[1 1 1]);

        set(h{k},'Marker','o','MarkerSize',MarkerSize,'MarkerFaceColor',series_color{k},'LineWidth',LineWidth);
        set(gca,'LineWidth',LineWidth)
        
        set(gca,'YTick',[-2:.1:2]);
%         set(gca,'YLim',[-.6 .5]);
        
        set(gca,'XTick',[0 7.5 15 22.5]);
        
        % set(gca,'XTick',[100 50 25 0]);
        % set(gca,'XLim',[-10 110]);
        set(gca,'XLim',[-2 17.5]);
        set(gca,'FontSize',FontSize)
    end
    legend([h{1}(1) h{2}(1) h{3}(1)],'opposite','same','blank','Location','NorthWest')
    legend boxoff;
end

box on;
xlabel('Time (s)','FontSize',24)
ylabel('fMRI response (% BOLD)','FontSize',24)
th = title([rois{1} 'No baseline subtraction. Threshold ='  num2str(cothresh) 'Phase window = ' num2str(phWindow(2)*(180/pi))]);
set(th,'FontAngle','italic');
set(th,'FontWeight','bold');
set(gcf,'Color','w');
set(gcf,'Position',[100 100 750 500]);
orient landscape
set(gcf,'PaperPosition',[0 0 11 8.5]);

% saveas(gcf,fullfile(figDir,['perceptualadaptation_',s{1},'_',num2str(cothresh*100),'_',rois{1}]),'fig');
saveas(gcf,fullfile(figDir,['perceptualadaptation_',s{1},'_',num2str(cothresh*100),'_',rois{1}]),'pdf');


%% plot results (No baseline subtraction - all scans separately)
% close all;
LineWidth = 1.5; %1.5;
MarkerSize = 8; %8; %10;
FontSize = 14; %14;
% timeList = 0:1.5:1.5*14; %14;

% if strcmp(subtype,'spscr')
%     series_color = {[0 0 1], [0 1 1]};
% elseif strcmp(subtype,'tmscr')
%     series_color = {[0 0 1], [1 0 0]};
% else
%     series_color = {[0 .5 0], [.5 0 0]};
% end

% series_color = {[1 1 0], [1 .5 0],[0 0 1],[1 1 0]}; % Same, different, blank
series_color = {[0 1 0], [1 0 0],[0 0 1],[1 1 0]};

series_style = {'-','--','--','-'};
figure, 
hold on
th = title([rois 'No baseline subtraction. Threshold ='  num2str(cothresh) 'Phase window = ' num2str(phWindow(2)*(180/pi))]);
set(th,'FontAngle','italic');
set(th,'FontWeight','bold');
set(gcf,'Color','w');
set(gcf,'Position',[100 100 750 500]);
orient landscape
set(gcf,'PaperPosition',[0 0 11 8.5]);

yh = ylabel('fMRI response (% BOLD)');
for j = 1:length(rois)
    for ns=1:nScans
    hold on, subplot(5,2,ns)
    for k = 1:nConds %[1 3] %these are conditions
        hold on
        eh = errorbar(timeList,singleCycleTs{k}{j}{ns},singleCycleStdErrTs{k}{j}{ns},series_style{k},'Color',[0 0 0],'LineWidth',LineWidth);
        plot([timeList(1) timeList(end)], [0 0], 'k--')

        h{k} = plot(timeList, singleCycleTs{k}{j}{ns},'o','MarkerEdgeColor',0*[1 1 1]);
        set(h{k},'Marker','o','MarkerSize',MarkerSize,'MarkerFaceColor',series_color{k},'LineWidth',LineWidth);
        set(gca,'LineWidth',LineWidth)

        set(gca,'YTick',[-2:.5:2]);
%         set(gca,'YLim',[-.6 .5]); 
        set(gca,'XTick',[0 7.5 15 22.5]);
        set(gca,'XLim',[-2 17.5]);
        set(gca,'FontSize',FontSize)
        box on;
    end
    end
    legend([h{1}(1) h{2}(1) h{3}(1)],'opposite','same','blank','Location','NorthWest')
    legend boxoff;
end


% xlabel('Time (s)','FontSize',24)
% ylabel('fMRI response (% BOLD)','FontSize',24)


% saveas(gcf,fullfile(figDir,['perceptualadaptation_allscans',s{1},'_',num2str(cothresh*100),'_',rois{1}]),'fig');
saveas(gcf,fullfile(figDir,['perceptualadaptation_allscans',s{1},'_',num2str(cothresh*100),'_',rois{1}]),'pdf');


%% plot results (Stimulus trial subtraction - old)
% 
% figDir = '/Volumes/Vision/MRI/PerceptualAdaptation/Figures';
% close all;
% 
% if nConds > 2
%     
%     LineWidth = 1.5;
%     MarkerSize = 8; %10;
%     FontSize = 14;
%     timeList = -10*1.5:1.5:1.5*14; %0:1.5:1.5*14; %14;
%     series_color = {[0 1 0], [1 0 0],[0 0 1]};
%     
%     series_style = {'-','--','---'};
%     
%     for j = 1:length(rois)
%         %     subplot(3,3,j) % no subplots for one ROI
%         for k = 1:2 %these are conditions
%             hold on
%             eh = errorbar(timeList,avgSingleCycleTs{k}{j}-aveAllStimCondsTs,avgSingleCycleStdErrTs{k}{j},series_style{k},'Color',[0 0 0],'LineWidth',LineWidth);
%             plot([timeList(1) timeList(end)], [0 0], 'k--')
%             %         removeErrorBarEnds(eh);
%             h{k} = plot(timeList, avgSingleCycleTs{k}{j}-aveAllStimCondsTs,'o','MarkerEdgeColor',0*[1 1 1]);
%             
%             set(h{k},'Marker','o','MarkerSize',MarkerSize,'MarkerFaceColor',series_color{k},'LineWidth',LineWidth);
%             set(gca,'LineWidth',LineWidth)
%             
%             set(gca,'YTick',[-.05:.1:.05]);
%             set(gca,'XTick',[0 7.5 15 22.5]);
%         set(gca,'YLim',[-.05 .05]);
%         set(gca,'FontSize',FontSize)
%             
%         end
%         legend([h{1}(1) h{2}(1)],'opposite','same','Location','SouthWest')
%         legend boxoff;
%     end
%     box on;
%     xlabel('Time (s)','FontSize',24)
%     ylabel('Relative fMRI response (% change in BOLD)','FontSize',24)
%     set(gcf,'Color','w');
%     set(gcf,'Position',[100 100 750 500]);
%     orient landscape
%     set(gcf,'PaperPosition',[0 0 11 8.5]);
%     th = title([rois 'Stimulus presence subtracted. Threshold ='  num2str(cothresh) 'Phase window = ' num2str(phWindow(2)*(180/pi))]);
%     set(th,'FontAngle','italic');
%     set(th,'FontWeight','bold');
%     saveas(gcf,fullfile(figDir,['perceptualadaptation_stimulussubtract_alltrialerrbars',s{1},num2str(cothresh*100),rois{1},num2str(numberOfUseableVoxels)]),'fig');
%     saveas(gcf,fullfile(figDir,['perceptualadaptation_stimulussubtract_alltrialerrbars',s{1},num2str(cothresh*100),rois{1},num2str(numberOfUseableVoxels)]),'pdf');
% end
% 

%% plot results (No baseline subtraction - Just average of all conditions - old)
% % figDir = '/Volumes/Vision/MRI/PerceptualAdaptation/Figures';
% % 
% % LineWidth = 1; %1.5;
% % MarkerSize = 5; %8; %10;
% % FontSize = 12; %14;
% % timeList = 0:1.5:1.5*14; %14;
% % 
% % series_color = {[0 0 0]};
% % 
% % series_style = {'-',};
% % hold off, figure, hold on
% % 
% %         eh = errorbar(timeList,aveAllCondsTs,stderrAllCondsTs,'-','Color',[0 0 0],'LineWidth',LineWidth);
% % 
% %         h{k} = plot(timeList, aveAllCondsTs,'o','MarkerEdgeColor',0*[1 1 1]);
% % 
% %         set(h{k},'Marker','o','MarkerSize',MarkerSize,'MarkerFaceColor',[0 0 0],'LineWidth',LineWidth);
% %         set(gca,'LineWidth',LineWidth)
% %         plot([4 4],[-0.35 0.3],'k--')
% %         plot([11.5 11.5],[-0.35 0.3],'k--')
% %         plot([19 19],[-0.35 0.3],'k--')
% %          
% %         set(gca,'YTick',[-1:.1:1]);
% %         set(gca,'YLim',[-.35 .3]);
% %         
% %         set(gca,'XTick',[0 4 7.5 11.5 15 19 22.5]);
% %         
% % 
% %         set(gca,'XLim',[-2 22.5]);
% %         set(gca,'FontSize',FontSize)
% % 
% % 
% % 
% % box on;
% % xlabel('Time (s)','FontSize',24)
% % ylabel('Average fMRI response (% BOLD)','FontSize',24)
% % th = title([rois 'Average of all conditions. Threshold ='  num2str(cothresh) 'Phase window = ' num2str(phWindow(2)*(180/pi))]);
% % set(th,'FontAngle','italic');
% % set(th,'FontWeight','bold');
% % set(gcf,'Color','w');
% % set(gcf,'Position',[100 100 750 500]);
% % orient landscape
% % set(gcf,'PaperPosition',[0 0 11 8.5]);
% % 
% % saveas(gcf,fullfile(figDir,['perceptualadaptationAverage__alltrialerrbars',s{1},'_',num2str(cothresh*100),'-',rois{1}]),'fig');
% % saveas(gcf,fullfile(figDir,['perceptualadaptationAverage__alltrialerrbars',s{1},'_',num2str(cothresh*100),'-',rois{1}]),'pdf');
% % 
% % 
% % 
% % 

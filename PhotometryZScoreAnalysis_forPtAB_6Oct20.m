%% Written by Khairunisa Ibrahim
% Adapted from TDT website
% For Moron-Concepcion's Lab WashU St Louis
% Analysis of photometry data to find zscore of photometry change
%% Clear workspace and close existing figures.
close all; clear all; clc;

%Add the SDK files to Matlab path
[MAINEXAMPLEPATH,name,ext] = fileparts(cd); % \TDTMatlabSDK\Examples
DATAPATH = fullfile(MAINEXAMPLEPATH, 'Examples'); % \TDTMatlabSDK\FolderName\
[SDKPATH,name,ext] = fileparts(MAINEXAMPLEPATH); % \TDTMatlabSDK
addpath(genpath(SDKPATH));

%% Set the directory of our data

BLOCKPATH = fullfile(DATAPATH,'N42-2-1-201006-115005'); %Data file have to be in the \TDTMatlabSDK\Examples\ directory

%To double check the blockpath, In Synapse, go to Menu > History. Find your block, then Right-Click > Copy path to clipboard

%% Setup the variables for the data you want to extract
% We will extract two different stream stores surrounding the 'Score' epoch event. 
% We are interested in a specific event code for the behavior onset.

REF_EPOC = 'PtAB'; %'Score'; % event store name. This holds behavioral codes that are read through Ports A & B on the front of the RZ
VALUE = 1; % score event code we are interested in
STREAM_STORE1 = 'x470A'; % name of the 470 store
STREAM_STORE2 = 'x405A'; % name of the 405 store
TRANGE = [-10 30]; % window size [start time relative to epoc onset, window duration]
BASELINE_PER = [-10 -1]; % baseline period within our window
ARTIFACT = Inf; % optionally set an artifact rejection level

% Now read the specified data from our block into a Matlab structure.
data = TDTbin2mat(BLOCKPATH, 'TYPE', {'epocs', 'scalars', 'streams'});

% Make a new epoc structure based on scores that is manually scored in OpenScope
% data.epocs.Score.onset = data.epocs.Cam1.notes.ts; 
% data.epocs.Score.offset = data.epocs.Score.onset + 0.05; 
% data.epocs.Score.name = 'Score';
% data.epocs.Score.data = double(data.epocs.Cam1.notes.index);

%% Use TDTfilter to extract data around our epoc event
% Using the 'TIME' parameter extracts data only from the time range around
% our epoc event. Use the 'VALUES' parameter to specify allowed values of
% the REF_EPOC to extract.  For stream events, the chunks of data are 
% stored in cell arrays structured as data.streams.(STREAM_STORE1).filtered

data = TDTfilter(data, REF_EPOC, 'VALUES', VALUE, 'TIME', TRANGE); %for PtAB

%% Optionally remove artifacts. If any waveform is above ARTIFACT level, or below -ARTIFACT level, remove it from the data set.
art1 = ~cellfun('isempty', cellfun(@(x) x(x>ARTIFACT), data.streams.(STREAM_STORE1).filtered, 'UniformOutput',false));
art2 = ~cellfun('isempty', cellfun(@(x) x(x<-ARTIFACT), data.streams.(STREAM_STORE1).filtered, 'UniformOutput',false));
good = ~art1 & ~art2;
data.streams.(STREAM_STORE1).filtered = data.streams.(STREAM_STORE1).filtered(good);

art1 = ~cellfun('isempty', cellfun(@(x) x(x>ARTIFACT), data.streams.(STREAM_STORE2).filtered, 'UniformOutput',false));
art2 = ~cellfun('isempty', cellfun(@(x) x(x<-ARTIFACT), data.streams.(STREAM_STORE2).filtered, 'UniformOutput',false));
good2 = ~art1 & ~art2;
data.streams.(STREAM_STORE2).filtered = data.streams.(STREAM_STORE2).filtered(good2);

numArtifacts = sum(~good) + sum(~good2);
%% Applying a time filter to a uniformly sampled signal means that the length of each segment could vary by one sample.  
% Let's find the minimum length so we can trim the excess off before calculating the mean.
minLength1 = min(cellfun('prodofsize', data.streams.(STREAM_STORE1).filtered));
minLength2 = min(cellfun('prodofsize', data.streams.(STREAM_STORE2).filtered));
data.streams.(STREAM_STORE1).filtered = cellfun(@(x) x(1:minLength1), data.streams.(STREAM_STORE1).filtered, 'UniformOutput',false);
data.streams.(STREAM_STORE2).filtered = cellfun(@(x) x(1:minLength2), data.streams.(STREAM_STORE2).filtered, 'UniformOutput',false);

% allSignals = cell2mat(data.streams.(STREAM_STORE2).filtered');
% This was from the example.
% allSignals470 = cell2mat(data.streams.(STREAM_STORE1).filtered');
% This was included in the downsampling section.
% allSignals405 = cell2mat(data.streams.(STREAM_STORE2).filtered');
% This was included in the downsampling section.

% downsample 1000x and average 470 signal
N = 100;
allSignals470 = cell2mat(data.streams.(STREAM_STORE1).filtered');
F470 = zeros(size(allSignals470(:,1:N:end-N+1)));
for ii = 1:size(allSignals470,1)
    F470(ii,:) = arrayfun(@(i) mean(allSignals470(ii,i:i+N-1)),1:N:length(allSignals470)-N+1);
end
minLength1 = size(F470,2);

% Create mean signal, standard error of signal, and DC offset of 470 signal
meanSignal470 = mean(F470);
stdSignal470 = std(double(F470))/sqrt(size(F470,1));
dcSignal470 = mean(meanSignal470);

% downsample 100x and average 405 signal
allSignals405 = cell2mat(data.streams.(STREAM_STORE2).filtered');
F405 = zeros(size(allSignals405(:,1:N:end-N+1)));
for ii = 1:size(allSignals405,1)
    F405(ii,:) = arrayfun(@(i) mean(allSignals405(ii,i:i+N-1)),1:N:length(allSignals405)-N+1);
end
minLength2 = size(F405,2);

% Create mean signal, standard error of signal, and DC offset of 405 signal
meanSignal405 = mean(F405);
stdSignal405 = std(double(F405))/sqrt(size(F405,1));
dcSignal405 = mean(meanSignal405);


%% Plot Epoch Averaged Response

% Create the time vector for each stream store
ts1 = TRANGE(1) + (1:minLength1) / data.streams.(STREAM_STORE1).fs*N;
ts2 = TRANGE(1) + (1:minLength2) / data.streams.(STREAM_STORE2).fs*N;

% Subtract DC offset to get signals on top of one another
meanSignal470 = meanSignal470 - dcSignal470;
meanSignal405 = meanSignal405 - dcSignal405;

% Plot the 470 and 405 average signals
figure;
subplot(2,2,1)
plot(ts1, meanSignal470, 'color',[0.4660, 0.6740, 0.1880], 'LineWidth', 3); hold on;
plot(ts2, meanSignal405, 'color',[0.8500, 0.3250, 0.0980], 'LineWidth', 3); 

% Plot vertical line at epoch onset, time = 0
line([0 0], [min(F470(:) - dcSignal470), max(F470(:)) - dcSignal470], 'Color', [.7 .7 .7], 'LineStyle','-', 'LineWidth', 3)

% Make a legend
legend('470 nm','405 nm','Onset', 'AutoUpdate', 'off');

% Create the standard error bands for the 405 signal
XX = [ts2, fliplr(ts2)];
YY = [meanSignal405 + stdSignal405, fliplr(meanSignal405 - stdSignal405)];

% Plot filled standard error bands.
h = fill(XX, YY, 'g');
set(h, 'facealpha',.25,'edgecolor','none')

% Repeat for 470
XX = [ts1, fliplr(ts1)];
YY = [meanSignal470 + stdSignal470, fliplr(meanSignal470 - stdSignal470)];
h = fill(XX, YY, 'r');
set(h, 'facealpha',.25,'edgecolor','none')

% Finish up the plot
axis tight
xlabel('Time, s','FontSize',12)
ylabel('V', 'FontSize', 12)
title(sprintf('Epoch Response, %d Trials (%d Artifacts Removed)', numel(data.streams.(STREAM_STORE1).filtered), numArtifacts))
set(gcf, 'Position',[100, 100, 800, 500])

% Heat Map based on z score of 405 fit subtracted 470

% Fitting 405 channel onto 470 channel to detrend signal bleaching
% Scale and fit data
% Algorithm sourced from Tom Davidson's Github:
% https://github.com/tjd2002/tjd-shared-code/blob/master/matlab/photometry/FP_normalize.m

bls = polyfit(F405(1:end), F470(1:end), 1);
Y_fit_all = bls(1) .* F405 + bls(2);
Y_dF_all = F470 - Y_fit_all;

zall = zeros(size(Y_dF_all));
for i = 1:size(Y_dF_all,1)
    ind = ts2(1,:) < BASELINE_PER(2) & ts2(1,:) > BASELINE_PER(1);
    zb = mean(Y_dF_all(i,ind)); % baseline period mean (-10sec to -6sec)
    zsd = std(Y_dF_all(i,ind)); % baseline period stdev
    zall(i,:)=(Y_dF_all(i,:) - zb)/zsd; % Z score per bin
end

% Standard error of the z-score
zerror = std(zall)/sqrt(size(zall,1));

%mean of the z-score
zmean = mean(zall);

% Plot heat map
subplot(2,2,2)

imagesc(ts2, 1, zall);
colormap('jet'); % c1 = colorbar; 
title(sprintf('Z-Score Heat Map, %d Trials (%d Artifacts Removed)', numel(data.streams.(STREAM_STORE1).filtered), numArtifacts));
ylabel('Trials', 'FontSize', 12);

% Fill band values for second subplot. Doing here to scale onset bar
% correctly
XX = [ts2, fliplr(ts2)];
YY = [mean(zall)-zerror, fliplr(mean(zall)+zerror)];

subplot(2,2,3)
plot(ts2, mean(zall), 'color',[0.8500, 0.3250, 0.0980], 'LineWidth', 3); hold on;
line([0 0], [min(YY), max(YY)], 'Color', [.7 .7 .7], 'LineWidth', 2)

h = fill(XX, YY, 'r');
set(h, 'facealpha',.25,'edgecolor','none')

% Finish up the plot
axis tight
xlabel('Time, s','FontSize',12)
ylabel('Z-score', 'FontSize', 12)
title(sprintf('470 nm Epoch Response, %d Trials (%d Artifacts Removed)', numel(data.streams.(STREAM_STORE1).filtered), numArtifacts))
%c2 = colorbar;

% Quantify changes as area under the curve for pre (-2 sec) and post (2 sec)
AUC=[]; % pre, post
AUC(1,1)=trapz(mean(zall(:,ts2(1,:) < 0 & ts2(1,:) > -2)));
AUC(1,2)=trapz(mean(zall(:,ts2(1,:) > 0 & ts2(1,:) < 2)));
subplot(2,2,4);
hBar = bar(AUC, 'FaceColor', [.8 .8 .8]);

% Run a two-sample T-Test
[h,p,ci,stats] = ttest2(mean(zall(:,ts2(1,:) < 0 & ts2(1,:) > -2)),mean(zall(:,ts2(1,:) > 0 & ts2(1,:) < 2)));

% Plot significance bar if p < .05
hold on;
centers = get(hBar, 'XData');
plot(centers(1:2), [1 1]*AUC(1,2)*1.1, '-k', 'LineWidth', 2)
p1 = plot(mean(centers(1:2)), AUC(1,2)*1.2, '*k');
set(gca,'xticklabel',{'Pre','Post'});
title({'Pre vs Post Response Changes', 'Area under curve'})
legend(p1, 'p < .05','Location','southeast');

set(gcf, 'Position',[100, 100, 900, 900])
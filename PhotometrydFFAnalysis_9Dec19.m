%% Written by Khairunisa Ibrahim
% Adapted from TDT website
% For Moron-Concepcion's Lab WashU St Louis
% Analysis of photometry data to find percentage change in delta F/F
%% Clear workspace and close existing figures.
close all; clear all; clc;

%Add the SDK files to Matlab path
[MAINEXAMPLEPATH,name,ext] = fileparts(cd); % \TDTMatlabSDK\Examples
DATAPATH = fullfile(MAINEXAMPLEPATH, 'Examples'); % \TDTMatlabSDK\FolderName\
[SDKPATH,name,ext] = fileparts(MAINEXAMPLEPATH); % \TDTMatlabSDK
addpath(genpath(SDKPATH));

%% Set the directory of our data

BLOCKPATH = fullfile(DATAPATH,'vGlut2-1-2-200605-141757'); %Data file have to be in the \TDTMatlabSDK\Examples\ directory

%To double check the blockpath, In Synapse, go to Menu > History. Find your block, then Right-Click > Copy path to clipboard

%% 
% Call the import function from the Matlab SDK
% https://www.tdt.com/support/sdk.html
data = TDTbin2mat(BLOCKPATH);

%%  Declare data stream and epoc names we will use downstream
% These are the field names for the relevant streams of the data struct
GCAMP = 'x470A';
ISOS = 'x405A';
%LICK = 'Ler_';
REF_EPOC = 'Score'; % event store name. This holds behavioral codes that are read through Ports A & B on the front of the RZ
VALUE = 1; % score event code we are interested in

% Make some pretty colors for later plotting
% http://math.loyola.edu/~loberbro/matlab/html/colorsInMatlab.html
black = [0.00,0.00,0.00];
red = [0.8500, 0.3250, 0.0980];
green = [0.36,0.63,0.01];
lightgreen = [0.58,0.80,0.42];
cyan = [0.3010, 0.7450, 0.9330];
gray1 = [.7 .7 .7];
gray2 = [.8 .8 .8];

%% Basic plotting and artifact removal
% Make a time array based on number of samples and sample freq of
% demodulated streams
time = (1:length(data.streams.(GCAMP).data))/data.streams.(GCAMP).fs;

%% Plot both unprocessed demodulated data streams
figure (1)
%'Position',[100, 100, 800, 400];
hold on;
p1 = plot(time, data.streams.(GCAMP).data,'color',green,'LineWidth',2);
p2 = plot(time, data.streams.(ISOS).data,'color',red,'LineWidth',2);
title('Raw Demodulated Responses','fontsize',16);
ylabel('mV','fontsize',16);
axis tight;
legend([p1 p2], {'GCaMP','ISOS'});

%% Artifact removal
% There is often a large artifact on the onset of LEDs turning on
% Remove data below a set time t
t = 8; % time threshold below which we will discard
ind = find(time>t,1); % find first index of when time crosses threshold
time = time(ind:end); % reformat vector to only include allowed time
data.streams.(GCAMP).data = data.streams.(GCAMP).data(ind:end);
data.streams.(ISOS).data = data.streams.(ISOS).data(ind:end);

%%  Plot at new time range without the artifact
figure (2)
%'Position',[100, 100, 800, 400];
hold on;
p1 = plot(time, data.streams.(GCAMP).data,'color',green,'LineWidth',2);
p2 = plot(time, data.streams.(ISOS).data,'color',red,'LineWidth',2);
title('Raw Demodulated Responses with Artifact Removed','fontsize',16);
xlabel('Seconds','fontsize',16)
ylabel('mV','fontsize',16);
axis tight;
legend([p1 p2], {'GCaMP','ISOS'});

%% Downsample data doing local averaging
% Average around every Nth point and downsample Nx

N = 1000; % multiplicative for downsampling
data.streams.(GCAMP).data = arrayfun(@(i)...
    mean(data.streams.(GCAMP).data(i:i+N-1)),...
    1:N:length(data.streams.(GCAMP).data)-N+1);
data.streams.(ISOS).data = arrayfun(@(i)...
    mean(data.streams.(ISOS).data(i:i+N-1)),...
    1:N:length(data.streams.(ISOS).data)-N+1);

%% Decimate time array and match length to demodulated stream
time = time(1:N:end);
time = time(1:length(data.streams.(GCAMP).data));

%% Detrending and dFF, same as Barker et al., 2017
bls = polyfit(data.streams.(ISOS).data,data.streams.(GCAMP).data,1);
Y_fit_all = bls(1) .* data.streams.(ISOS).data + bls(2);
Y_dF_all = data.streams.(GCAMP).data - Y_fit_all; %dF (units mV) is not dFF

%% Full dFF according to Lerner et al. 2015
% http://dx.doi.org/10.1016/j.cell.2015.07.014
% dFF using 405 fit as baseline
dFF = 100*(Y_dF_all)./Y_fit_all; %dividing each element by the corresponding element in ISOS
std_dFF = std(double(dFF));
%% To Find z score of dFF
dFF_raw = (Y_dF_all)./Y_fit_all;
norm_dFF = zscore(dFF_raw);

%% Plot dFF 
figure (3)
%'Position',[100, 100, 800, 400];
p1 = plot(time, dFF, 'Color',green,'LineWidth',2); hold on;
%p2 = plot(LICK_x, y_scale*(LICK_y) + y_shift,'color',cyan,'LineWidth',2);
title('Detrended, y-shifted dFF','fontsize',16);
legend([p1],'GCaMP')
%legend([p1 p2],'GCaMP','Lick Epoc');
axis tight
%% Finding peaks
%Find the peaks are at least 2 dFF.
[peakswthreshold] = findpeaks(dFF,time,'MinPeakHeight',2);
[meanPeakswThreshold] = mean(peakswthreshold);

%Find the peaks that drop at least 1.5 dFF on either side before the signal attains a higher value.
[peaks] = findpeaks(dFF,time,'MinPeakProminence',1.5);
[meanPeaks] = mean(peaks);

%Find the peaks that are at least 1.5 dFF higher than their neighboring samples.
[peakswNeigh] = findpeaks(dFF,time,'Threshold',1.5);

%Find the peaks that are at least 1.5 dFF higher than their neighboring
%samples with width reference at half height
[peakswWidth] = findpeaks(dFF,time,'Threshold',1.5,'Annotate','extents','WidthReference','halfheight');
[meanPeakswWidth] = mean(peakswWidth);

% findpeaks that are at least 2 dFF and have a width of at least 0.5 sec at
% half amplitude
%[peakswAmpandWidth] = (dFF,time,'MinPeakProminence',1.5,'WidthReference','halfprom','MinPeakWidth',0.5);
% title('Prom of 1.5, Width of 0.5','fontsize',16);


%% make a new epoc structure based on scores
data.epocs.Score.onset = data.epocs.Cam1.notes.ts; 
data.epocs.Score.offset = data.epocs.Score.onset + 0.05; 
data.epocs.Score.name = 'Score';
data.epocs.Score.data = double(data.epocs.Cam1.notes.index);

%% %% Use TDTfilter to extract data around our epoc event
% Using the 'TIME' parameter extracts data only from the time range around
% our epoc event. Use the 'VALUES' parameter to specify allowed values of
% the REF_EPOC to extract.  For stream events, the chunks of data are 
% stored in cell arrays structured as data.streams.(STREAM_STORE1).filtered
data = TDTfilter(data, REF_EPOC, 'VALUES', VALUE); %, 'TIME', TRANGE);

%% Time Filter around the Epochs
%Note that this is using dFF of the full time-series and not peri-event dFF where f0 is taken from a pre-event baseline period. 

PRE_TIME = 15; % 15 seconds before event onset
POST_TIME = 30; % 30 seconds after
fs = data.streams.(GCAMP).fs/N; % recall we downsampled by N = 1000 earlier
% time span for peri-event filtering, PRE and POST
TRANGE = [-1*PRE_TIME*floor(fs),POST_TIME*floor(fs)];

%% creating empty cell or matrix to store the values
trials = numel(data.epocs.Score.data);
dFF_snips = cell(trials,1);
array_ind = zeros(trials,1);
pre_stim = zeros(trials,1);
post_stim = zeros(trials,1);

%% Trim the stream based on the trigger onset
for i = 1:trials
    % If the bout cannot include pre-time seconds before event, make zero
    if data.epocs.Score.onset(i) < PRE_TIME
        dFF_snips{i} = single(zeros(1,(TRANGE(2)-TRANGE(1))));
        continue
    else
        % Find first time index after bout onset
        array_ind(i) = find(time > data.epocs.Score.onset(i),1);
       
        % Find index corresponding to pre and post stim durations
        pre_stim(i) = array_ind(i) + TRANGE(1);
        post_stim(i) = array_ind(i) + TRANGE(2);
        dFF_snips{i} = dFF(pre_stim(i):post_stim(i));
    end
end

%% Trim the snippets the same cells to the min length
minLength = min(cellfun('prodofsize', dFF_snips));
dFF_snips = cellfun(@(x) x(1:minLength), dFF_snips, 'UniformOutput',false);

% Convert to a matrix and get mean
allSignals = cell2mat(dFF_snips);
mean_allSignals = mean(allSignals);
std_allSignals = std(mean_allSignals);

% Make a time vector snippet for peri-events
peri_time = (1:length(mean_allSignals))/fs - PRE_TIME;

%% Plot the Per-event stimulus
% Make a standard deviation fill for mean signal
figure (4)
%('Position',[100, 100, 600, 750])
%subplot(2,1,1)
xx = [peri_time, fliplr(peri_time)];
yy = [mean_allSignals + std_allSignals,...
    fliplr(mean_allSignals - std_allSignals)];
h = fill(xx, yy, 'g'); % plot this first for overlay purposes
hold on;
set(h, 'facealpha', 0.25, 'edgecolor', 'none');

% Set specs for min and max value of event line.
% Min and max of either std or one of the signal snip traces
linemin = min(min(min(allSignals)),min(yy));
linemax = max(max(max(allSignals)),max(yy));

% Plot the line next
l1 = line([0 0], [linemin, linemax],...
    'color','black', 'LineStyle', '-', 'LineWidth', 2);
% Plot the signals and the mean signal
p1 = plot(peri_time, allSignals', 'color', gray1);
p2 = plot(peri_time, mean_allSignals, 'color', green, 'LineWidth', 3);
hold off;

% Make a legend and do other plot things
legend([l1, p1(1), p2, h],...
    {'Onset','Trial Traces','Mean Response','Std'},...
    'Location','northeast');
title('Peri-Event Trial Responses','fontsize',16);
ylabel('\DeltaF/F','fontsize',16);
grid
axis tight;
% Make an invisible colorbar so this plot aligns with one below it
temp_cb = colorbar('Visible', 'off');

%% Plot the heat map
% Heat map
figure (5)
%subplot(2,1,2)
imagesc(peri_time, 1, allSignals); % this is the heatmap
set(gca,'YDir','normal') % put the trial numbers in better order on y-axis
colormap(gray()) % colormap otherwise defaults to perula
title('Bout Heat Map','fontsize',16)
ylabel('Trial Number','fontsize',16)
xlabel('Seconds from onset','fontsize',16)
cb = colorbar;
ylabel(cb, 'dFF','fontsize',16)
axis tight;

%% Calculating Area Under the Curve (AUC)

%taking from the y-axis
AUC = [];
AUC(1,1)=trapz(mean(allSignals(:,peri_time(1,:) < 0 & peri_time(1,:) > -5)));
AUC(1,2)=trapz(mean(allSignals(:,peri_time(1,:) > 0 & peri_time(1,:) < 5)));
%subplot(4,1,4);
figure (6)
hBar = bar(AUC, 'FaceColor', [.8 .8 .8]);
%% Calculating Area Under the Curve (AUC) absolute value

%AUCabs = [];
%AUCabs(1,1)=trapz(abs(mean_allSignals(:,peri_time(1,:) < 0 & peri_time(1,:) > -5)));
%AUCabs(1,2)=trapz(abs(mean_allSignals(:,peri_time(1,:) > 0 & peri_time(1,:) < 5)));
%subplot(4,1,4);
%figure (7)
%hBar = bar(AUCabs, 'FaceColor', [.8 .8 .8]);
%% Plotting normalized dF/F 

figure (7)
plot(time,norm_dFF)
ylabel('Normalized \Delta F/F')
xlabel('Time (Seconds)')
title('Normalized \Delta F/F for Recording ')

figure (8)
xx = [peri_time, fliplr(peri_time)];
yy = [mean_allSignals + std_allSignals,...
    fliplr(mean_allSignals - std_allSignals)];
h = fill(xx, yy, 'g'); % plot this first for overlay purposes
hold on;
set(h, 'facealpha', 0.25, 'edgecolor', 'none');

% Set specs for min and max value of event line.
% Min and max of either std or one of the signal snip traces
linemin = min(min(min(allSignals)),min(yy));
linemax = max(max(max(allSignals)),max(yy));

% Plot the line next
l1 = line([0 0], [linemin, linemax],...
    'color','black', 'LineStyle', '-', 'LineWidth', 2);
% Plot the signals and the mean signal
p1 = plot(peri_time, allSignals', 'color', gray1);
p2 = plot(peri_time, mean_allSignals, 'color', green, 'LineWidth', 3);
hold off;

% Make a legend and do other plot things
legend([l1, p1(1), p2, h],...
    {'Onset','Trial Traces','Mean Response','Std'},...
    'Location','northeast');
title('Peri-Event Trial Responses','fontsize',16);
ylabel('\DeltaF/F','fontsize',16);
grid
axis tight;
% Make an invisible colorbar so this plot aligns with one below it
temp_cb = colorbar('Visible', 'off');




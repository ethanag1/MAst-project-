clear 
clc
%maxNumCompThreads(feature('numcores'))

%% setup
dir_coi = '/home/ethan/coi'; % set filepath to coi toolbox on pc
cd(dir_coi) % change matlab's current directory to it
coi_init % make coi code available for matlab to use
% make sure to also add the gcmi toolbox to path
addpath('/home/ethan/gcmi/matlab');
savepath
addpath('/home/ethan/scripts_report'); % Just add the folder
savepath
% loads eeglab 
run('/home/ethan/scripts/eeglab2021.0/eeglab.m')

% load the data 
load('/home/ethan/Data/BOD_volitional_streaming_04.mat')

%% finding the indices where perception is 1 stream or 2 streams

% Find indices where perception is 1 stream
oneIndices = find(perception == 1);

% Find indices where perception is 2 streams
twoIndices = find(perception == 2);

% Display results
disp(['Length of oneIndices: ', num2str(length(oneIndices))]);
disp(['Length of twoIndices: ', num2str(length(twoIndices))]);

%% find the corresponding values in stim for these indices in oneIndices and twoIndices

% Initialize vectors with 3s, same length as stim
stimOne = 3 * ones(size(stim)); 
stimTwo = 3 * ones(size(stim));

% found values in stim retain the same indices
stimOne(oneIndices) = stim(oneIndices); % Retain original indices
stimTwo(twoIndices) = stim(twoIndices); % Retain original indices

%% identifying the indices at which the start of each tone occurs in stimOne and stimTwo (i.e. the start of each block of 1s in stimOne and stimTwo)

% Find indices where stimOne is 1
oneIndices_stimOne = find(stimOne == 1);

% Identify breaks between separate sections
diffs = diff(oneIndices_stimOne);
breaks = find(diffs > 1);

% Ensure startIndicesOne is defined properly
if isempty(breaks)  % Case where all 1s are in one continuous block
    startIndices_stimOne = oneIndices_stimOne(1);
else
    startIndices_stimOne = [oneIndices_stimOne(1); oneIndices_stimOne(breaks + 1)];  
end

% Ensure that startIndices is a column vector
startIndices_stimOne = startIndices_stimOne(:);

% Find indices where stimTwo is 1
oneIndices_stimTwo = find(stimTwo == 1);

% Identify breaks between separate sections
diffs = diff(oneIndices_stimTwo);
breaks = find(diffs > 1);

% Ensure startIndices_stimTwo is defined properly
if isempty(breaks)  % Case where all 1s are in one continuous block
    startIndices_stimTwo = oneIndices_stimTwo(1);
else
    startIndices_stimTwo = [oneIndices_stimTwo(1); oneIndices_stimTwo(breaks + 1)];  
end

% Ensure that startIndices is a column vector
startIndices_stimTwo = startIndices_stimTwo(:);

disp('1 stream tones');
disp(length(startIndices_stimOne));
disp('2 streams tones');
disp(length(startIndices_stimTwo));

%% Identifying the indices of the start of each triplet
% i.e. since each triplet contains 3 tones need to take every 3rd value in
% the vector "startIndices" starting from the first one

% Select every third index starting from the first one
selectedIndices_stimOne = startIndices_stimOne(1:3:end-2); % note: subtracted 2 from the end otherwise it will include/count 
selectedIndices_stimTwo = startIndices_stimTwo(1:3:end-2); % only part a triplet so subtracting 2 ensures that only full triplets are counted

disp('1 stream triplets');
disp(length(selectedIndices_stimOne));
disp('2 streams triplets');
disp(length(selectedIndices_stimTwo));
%% Filter data using EEGlab

% Channels on first dimension
dp_data = dp_data'; 

% import data to EEGlab format
EEG = pop_importdata('dataformat','matlab','nbchan',0,'data',dp_data,'setname','data','srate',1200,'pnts',0,'xmin',0);

% Filter data 0.5-40 Hz
EEG = pop_eegfiltnew(EEG, 'locutoff',0.5,'hicutoff',40,'plotfreqz',1);

% back to channels second dimension
dp_data = EEG.data';

%% creating sections of length -60 entries (50ms before onset of triplet) and +660 entries 
% (550ms i.e. the length of the triplet) for each entry in "selectedIndices" and making cuts 
% in dp_data and transforming the 2D dp_data matrix into a 3d one with the
% 3rd dim being each section cut from the dp_data matrix

% Define the range around each selected index
preTime = 60;    % Number of points before the selected index (baseline period)
postTime = 660;  % Number of points after the selected index (triplet period)
sectionLength = preTime + postTime + 1; % Total length of each section

% Define the number of columns in the dp_data 2D matrix
numChannels = size(dp_data, 2);
% define the number of rows in dp_data 2D matrix 
numTimePoints = size(dp_data,1);

% Number of sections based on selected indices
numSections1 = length(selectedIndices_stimOne);

% Initialize 3D matrix to store extracted sections
BaselineMatrix_3D_1stream = zeros(sectionLength, numChannels, numSections1);

% Extract sections
for i = 1:numSections1
    % Define start and end indices for each section
    centerIdx = selectedIndices_stimOne(i);
    startIdx = centerIdx - preTime;
    endIdx = centerIdx + postTime;

    % Ensure the indices stay within bounds
    if startIdx >= 1 && endIdx <= numTimePoints
        BaselineMatrix_3D_1stream(:, :, i) = dp_data(startIdx:endIdx, :); % this might be wrong, need to check this 
    else
        warning('Skipping section %d: Indices out of bounds.', i);
    end
end

% Number of sections based on selected indices
numSections2 = length(selectedIndices_stimTwo);

% Initialize 3D matrix to store extracted sections
BaselineMatrix_3D_2streams = zeros(sectionLength, numChannels, numSections2);

% Extract sections
for i = 1:numSections2
    % Define start and end indices for each section
    centerIdx = selectedIndices_stimTwo(i);
    startIdx = centerIdx - preTime;
    endIdx = centerIdx + postTime;

    % Ensure the indices stay within bounds
    if startIdx >= 1 && endIdx <= numTimePoints
        BaselineMatrix_3D_2streams(:, :, i) = dp_data(startIdx:endIdx, :); % this might be wrong, need to check this 
    else
        warning('Skipping section %d: Indices out of bounds.', i);
    end
end

% Display confirmation
disp('✅ Successfully transformed the 2D matrix into a 3D matrix.');
disp(['Size of 3D matrix 1 stream: ', num2str(size(BaselineMatrix_3D_1stream,1)), ' x ', num2str(size(BaselineMatrix_3D_1stream,2)), ' x ', num2str(size(BaselineMatrix_3D_1stream,3))]);
disp(['Size of 3D matrix 2 streams: ', num2str(size(BaselineMatrix_3D_2streams,1)), ' x ', num2str(size(BaselineMatrix_3D_2streams,2)), ' x ', num2str(size(BaselineMatrix_3D_2streams,3))]);

%% Calculating the mean value of the first 60 entries/rows in BaselineMatrix_3D
% (i.e. the mean baseline value) for each 3rd dim entry separately and
% subtracting that from every entry for the corresponding 3rd dim entry for
% each channel separately 

% Compute mean over the first 60 rows for each channel and section (the
% baseline)
baselineMean_1stream = mean(BaselineMatrix_3D_1stream(1:60, :, :), 1); % Size: [1, numChannels, numSections]
baselineMean_2streams = mean(BaselineMatrix_3D_2streams(1:60, :, :), 1); % Size: [1, numChannels, numSections]
% Subtract the baseline mean from all rows
Triplet_baselinecorrected_Matrix_1stream = BaselineMatrix_3D_1stream - baselineMean_1stream;  % Broadcasting in MATLAB
Triplet_baselinecorrected_Matrix_2streams = BaselineMatrix_3D_2streams - baselineMean_2streams;  % Broadcasting in MATLAB
disp('✅ Successfully subtracted baselinemean from each colomn and section.');

%% finding the indices where perception is 1 stream or 2 streams

% Find indices where perception is 1 stream
oneIndices = find(perception == 1);

% Find indices where perception is 2 streams
twoIndices = find(perception == 2);

% Display results
disp(['Length of oneIndices: ', num2str(length(oneIndices))]);
disp(['Length of twoIndices: ', num2str(length(twoIndices))]);

%% find the corresponding values in stim for these indices in oneIndices and twoIndices

% Initialize vectors with 3s, same length as stim
stimOne = 3 * ones(size(stim)); 
stimTwo = 3 * ones(size(stim));

% found values in stim retain the same indices
stimOne(oneIndices) = stim(oneIndices); % Retain original indices
stimTwo(twoIndices) = stim(twoIndices); % Retain original indices

%% identifying the indices at which the start of each tone occurs in stimOne and stimTwo (i.e. the start of each block of 1s in stimOne and stimTwo)

% Find indices where stimOne is 1
oneIndices_stimOne = find(stimOne == 1);

% Identify breaks between separate sections
diffs = diff(oneIndices_stimOne);
breaks = find(diffs > 1);

% Ensure startIndicesOne is defined properly
if isempty(breaks)  % Case where all 1s are in one continuous block
    startIndices_stimOne = oneIndices_stimOne(1);
else
    startIndices_stimOne = [oneIndices_stimOne(1); oneIndices_stimOne(breaks + 1)];  
end

% Ensure that startIndices is a column vector
startIndices_stimOne = startIndices_stimOne(:);

% Find indices where stimTwo is 1
oneIndices_stimTwo = find(stimTwo == 1);

% Identify breaks between separate sections
diffs = diff(oneIndices_stimTwo);
breaks = find(diffs > 1);

% Ensure startIndices_stimTwo is defined properly
if isempty(breaks)  % Case where all 1s are in one continuous block
    startIndices_stimTwo = oneIndices_stimTwo(1);
else
    startIndices_stimTwo = [oneIndices_stimTwo(1); oneIndices_stimTwo(breaks + 1)];  
end

% Ensure that startIndices is a column vector
startIndices_stimTwo = startIndices_stimTwo(:);

disp('1 stream tones');
disp(length(startIndices_stimOne));
disp('2 streams tones');
disp(length(startIndices_stimTwo));

%% Identifying the indices of the start of each triplet
% i.e. since each triplet contains 3 tones need to take every 3rd value in
% the vector "startIndices" starting from the first one

% Select every third index starting from the first one
selectedIndices_stimOne = startIndices_stimOne(1:3:end-2); % note: subtracted 2 from the end otherwise it will include/count 
selectedIndices_stimTwo = startIndices_stimTwo(1:3:end-2); % only part a triplet so subtracting 2 ensures that only full triplets are counted

disp('1 stream triplets');
disp(length(selectedIndices_stimOne));
disp('2 streams triplets');
disp(length(selectedIndices_stimTwo));

%% combining the two types of triplets identified into 1 single vector retaining the correct trial ordering

% Combine and sort the indices
combinedIndices = sort([selectedIndices_stimOne; selectedIndices_stimTwo]);

%rename the combined indices to selectedIndices
selectedIndices = combinedIndices;

%% Filter data using EEGlab

% Channels on first dimension
dp_data = dp_data'; 

% import data to EEGlab format
EEG = pop_importdata('dataformat','matlab','nbchan',0,'data',dp_data,'setname','data','srate',1200,'pnts',0,'xmin',0);

% Filter data 0.5-40 Hz
EEG = pop_eegfiltnew(EEG, 'locutoff',0.5,'hicutoff',40,'plotfreqz',1);

% back to channels second dimension
dp_data = EEG.data';

%% creating sections of length +660 entries 
% (550ms i.e. the length of the triplet) for each entry in "selectedIndices" and making cuts 
% in dp_data and transforming the 2D dp_data matrix into a 3d one with the
% 3rd dim being each section cut from the dp_data matrix

% Define the range around each selected index
preTime = 60;    % Number of points before the selected index (baseline period)
postTime = 660;  % Number of points after the selected index (triplet period)
sectionLength = preTime + postTime + 1; % Total length of each section

% Define the number of columns in the dp_data 2D matrix
numChannels = size(dp_data, 2);
% define the number of rows in dp_data 2D matrix 
numTimePoints = size(dp_data,1);

% Number of sections based on selected indices
numSections = length(selectedIndices);

% Initialize 3D matrix to store extracted sections
BaselineMatrix_3D = zeros(sectionLength, numChannels, numSections);

% Extract sections
for i = 1:numSections
    % Define start and end indices for each section
    centerIdx = selectedIndices(i);
    startIdx = centerIdx - preTime;
    endIdx = centerIdx + postTime;

    % Ensure the indices stay within bounds
    if startIdx >= 1 && endIdx <= numTimePoints
        BaselineMatrix_3D(:, :, i) = dp_data(startIdx:endIdx, :); % this might be wrong, need to check this 
    else
        warning('Skipping section %d: Indices out of bounds.', i);
    end
end

% Display confirmation
disp('✅ Successfully transformed the 2D matrix into a 3D matrix.');
disp(['Size of 3D matrix: ', num2str(size(BaselineMatrix_3D,1)), ' x ', num2str(size(BaselineMatrix_3D,2)), ' x ', num2str(size(BaselineMatrix_3D,3))]);

%% Calculating the mean value of the first 60 entries/rows in BaselineMatrix_3D
% (i.e. the mean baseline value) for each 3rd dim entry separately and
% subtracting that from every entry for the corresponding 3rd dim entry for
% each channel separately 

% Compute mean over the first 60 rows for each channel and section (the
% baseline)
baselineMean = mean(BaselineMatrix_3D(1:60, :, :), 1); % Size: [1, numChannels, numSections]

% Subtract the baseline mean from all rows
Triplet_baselinecorrected_Matrix = BaselineMatrix_3D -  baselineMean;   

disp('✅ Successfully subtracted baselinemean from each colomn and section.');

%% trimming down the matrix to only include the triplet period (i.e. the last 550 entries/rows in the matrix)
% and not the baseline period

% Exclude the first 60 timepoints
BaselineMatrix_3D_trimmed = Triplet_baselinecorrected_Matrix(:, :, :);

% without baseline correction
%BaselineMatrix_3D_trimmed = BaselineMatrix_3D(61:end, :, :);
% Display the new size to confirm
disp('New size of the trimmed matrix:');
disp(size(BaselineMatrix_3D_trimmed));

%% permuting the inputs so that the matrix is trials x samples x channels

BaselineMatrix_3D_ep = permute(BaselineMatrix_3D_trimmed, [3 1 2]);

% Display confirmation
disp('✅ Successfully permuted matrix.');
disp(['Size of 3D matrix: ', num2str(size(BaselineMatrix_3D_ep,1)), ' x ', num2str(size(BaselineMatrix_3D_ep,2)), ' x ', num2str(size(BaselineMatrix_3D_ep,3))]);

%% testing parameters
% edit these to change:
    % how many timepoints you compute mi/coi on
    % how many permutations you run for significance testing
    n_t = 721; % number of timepoints
    n_pl = 100; % number of permutations
%% creating vector that labels the trials in the input according to the stimulus

% Initialize l_ vector with zeros
l_ = zeros(size(combinedIndices));

% Label trials: 0 for stimOne, 1 for stimTwo
l_(ismember(combinedIndices, selectedIndices_stimOne)) = 0;
l_(ismember(combinedIndices, selectedIndices_stimTwo)) = 1;
disp('length of l_ vector created');
disp(length(l_));
disp(l_);
%% extract inputs

% to neatly plot the mi, we also need the vector labelling timepoints
t_ = -60:1:660;
n_t = min(n_t, numel(t_)); % this stops selection of more timepoints than there are in the input
t_ = t_(1:n_t)/1.2; % convert to seconds

%%
% we can decode the mi for each timepoint, and its statistical significance, 
% by inputting those arguments to mi_decode 

% Get the number of channels and timepoints
[num_trials, num_timepoints, num_channels] = size(BaselineMatrix_3D_ep);

% Initialize storage for results as numeric matrices
mi_all = zeros(num_timepoints, num_channels);
zmi_all = zeros(num_timepoints, num_channels);

% Loop over each channel
for ch = 1:num_channels
    % Extract data for the current channel
    rr_ = BaselineMatrix_3D_ep(:, :, ch);
    
    % Compute MI and statistical significance
    [mi_, zmi_] = mi_decode(l_, rr_, 't', n_t, 'pl', n_pl);
    
    % Store results in pre-allocated matrices
    mi_all(:, ch) = mi_;
    zmi_all(:, ch) = zmi_;
end

% Output confirmation
disp('Mutual information computed for all channels.');

%% Identify channels with significant MI values

% Threshold for statistical significance (1.645 corresponds to p < 0.05, one-tailed)
threshold = 1.645;

% Find channels where any time point has significant zMI
significant_channels = find(any(zmi_all > threshold, 1));
%significant_channels = [10, 22, 23, 24, 25, 41, 44, 45, 46, 47, 48, 76, 77, 93, 149, 155, 160, 161, 169, 170, 171];

% Display results

disp('Significant channels:');
disp(significant_channels);
disp('number of channels with significant MI values');
disp(length(significant_channels));


%% Plot of single channel with error bars and MI underneath

% Choose which column to plot 
columnIndex = 92;

% Extract the selected column from both matrices
columnData_1stream = Triplet_baselinecorrected_Matrix_1stream(:, columnIndex);
columnData_2streams = Triplet_baselinecorrected_Matrix_2streams(:, columnIndex);

% Channels in first dimension for both matrices
dp_data_ep_1stream = permute(Triplet_baselinecorrected_Matrix_1stream, [2 1 3]);
dp_data_ep_2streams = permute(Triplet_baselinecorrected_Matrix_2streams, [2 1 3]);

% Create the x-axis 
x = -60:length(columnData_1stream)-61;
time = x / srate * 1000;

% Define x-axis ticks
xTicks = [0 250 500];

% Compute symmetric y-axis limits
yLimits = [-25, 25];

% Compute x-axis limits
xLimits = [min(time), max(time)];

% Define timepoints in milliseconds for MI
channel_to_plot = columnIndex; % Channel to plot
num_timepoints = size(mi_all, 1);
sampling_rate = 1200;
time_ms = ((1:num_timepoints) - 60) / sampling_rate * 1000; % Shift by -60 samples

% Get full MI and ZMI values
mi_channel = mi_all(:, channel_to_plot);
zmi_channel = zmi_all(:, channel_to_plot);
max_mi_val = max(mi_channel);
min_mi_val = min(mi_channel); % <- new minimum
mi_yticks = [0, max_mi_val]; %max_mi_val

% Create tiled layout
fig = figure;
t = tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
fig.Position(3:4) = [600 600];  % Square figure

%% ERP Plot
nexttile;
hold on;
stdshade_acj(squeeze(dp_data_ep_1stream(columnIndex,:,:))', 0.2, 'r', time, 2.5);
stdshade_acj(squeeze(dp_data_ep_2streams(columnIndex,:,:))', 0.2, 'b', time, 2.5);
hold off;

ylabel('Voltage (\muV) ', 'FontSize', 16, 'FontWeight', 'bold');

ax1 = gca;
ax1.YLim = yLimits;
ax1.XLim = xLimits;
ax1.XTick = xTicks;
ax1.XTickLabel = {};  % Remove x-tick labels but keep ticks
ax1.YTick = [-25 0 25];
ax1.XAxisLocation = 'origin';
ax1.YAxisLocation = 'left';
ax1.Box = 'off';
ax1.TickDir = 'out';
ax1.LineWidth = 2;  % Bolder axis
ax1.FontSize = 20;
ax1.TickLength = [0.015 0.015];

% Vertical dashed line at x = 0
hold on;
h = plot([0 0], yLimits, '--k', 'LineWidth', 1);
h.Color = [0, 0, 0, 0.3];

% Legend
%h1 = plot(NaN, NaN, 'r', 'LineWidth', 2.5);
%h2 = plot(NaN, NaN, 'b', 'LineWidth', 2.5);
%legend([h1, h2], {'1-stream', '2-streams'}, 'Location', 'northeast', 'FontSize', 13);
hold off;

%% MI Plot
nexttile;
hold on;
plot(time_ms, mi_channel, 'k', 'LineWidth', 2.5);

% Significant time points
significant_points = find(zmi_channel > 1.645);
if ~isempty(significant_points)
    y_max = max_mi_val * 1.07;
    for t = significant_points'
        plot([time_ms(t) time_ms(t)], [y_max, y_max * 1.07], 'Color', [0.5 0 0.5], 'LineWidth', 1.1);
    end
end

% Vertical dashed line at x = 0
h = plot([0 0], [min_mi_val, max_mi_val*1.1], '--k', 'LineWidth', 1);
h.Color = [0, 0, 0, 0.3];

xlabel('Time (ms)', 'FontSize', 16, 'FontWeight', 'bold');
ylabel('MI (bits)', 'FontSize', 16, 'FontWeight', 'bold');

ax2 = gca;
ax2.XLim = xLimits;
ax2.XTick = xTicks;
ax2.YTick = mi_yticks;
ax2.YLim = [min_mi_val, max_mi_val * 1.145]; % <- Adjust y-limits to remove empty space
ax2.LineWidth = 2;  % Bolder axis
ax2.Box = 'off';
ax2.FontSize = 20;
ax2.TickDir = 'out';
ax2.TickLength = [0.015 0.015];
hold off;

% Save the combined figure
saveas(fig, '/home/ethan/fig_report/ERP_separated_triplets_with_MI_plot_uncorrected_for_channel92_BOD(p2).png');


%columnIndex = 92;
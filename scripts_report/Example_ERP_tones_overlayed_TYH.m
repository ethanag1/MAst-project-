% loads eeglab 
run('/home/ethan/scripts/eeglab2021.0/eeglab.m')

% load the data
load('/home/ethan/Data/TYH_volitional_streaming_04.mat')
whos

%% identifying the indices at which the start of each tone occurs in stim (i.e. the start of each block of 1s in stim)

% Find indices where stim is 1
oneIndices = find(stim == 1);

% Identify breaks between separate sections
diffs = diff(oneIndices);
breaks = find(diffs > 1);

% Ensure startIndices is defined properly
if isempty(breaks)  % Case where all 1s are in one continuous block
    startIndices = oneIndices(1);
else
    startIndices = [oneIndices(1); oneIndices(breaks + 1)];  
end

% Ensure that startIndices is a column vector
startIndices = startIndices(:);

%disp(startIndices);
disp(length(startIndices));

%% Identifying the indices of the start of each triplet
% i.e. since each triplet contains 3 tones need to take every 3rd value in
% the vector "startIndices" starting from the first one

% Select every third index starting from the first one
selectedIndices = startIndices(1:3:end);

% Display the selected indices
%disp('Selected indices:');
%disp(selectedIndices);
% "selectedIndices" should have a length 500, since "startIndices" has
% length 1500
disp(length(selectedIndices));

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
Triplet_baselinecorrected_Matrix = BaselineMatrix_3D - baselineMean;  % Broadcasting in MATLAB

disp('✅ Successfully subtracted baselinemean from each colomn and section.');

%% averaging over the 3rd dim to obtain a 2D matrix again

% Compute the mean over the 3rd dimension (averaging across sections)
Triplet_baselinecorrected_Matrix_avrg = mean(Triplet_baselinecorrected_Matrix, 3);

% Ensure the result is truly 2D by removing the singleton dimension
Triplet_baselinecorrected_Matrix_avrg = squeeze(Triplet_baselinecorrected_Matrix_avrg);

% Display confirmation
disp('✅ Successfully collapsed the 3D matrix into a true 2D matrix by averaging.');
disp(['Size of new 2D matrix: ', num2str(size(Triplet_baselinecorrected_Matrix_avrg,1)), ' x ', num2str(size(Triplet_baselinecorrected_Matrix_avrg,2))]);


%% Plot of single channel with error bars (cropped to post-stimulus only, annotated with stacked tone bars)

% Choose which column to plot 
columnIndex = 171;

% Extract post-stimulus portion only
columnData = Triplet_baselinecorrected_Matrix_avrg(61:end, columnIndex);

% Channels in first dimension
dp_data_ep = permute(Triplet_baselinecorrected_Matrix, [2 1 3]);

% Extract corresponding post-stimulus data
dp_data_ep_post = dp_data_ep(columnIndex, 61:end, :);

% Create time vector from 0 ms onward
x = 0:size(columnData, 1) - 1;
time = x / srate * 1000;

% Plot
figure;
stdshade_acj(squeeze(dp_data_ep_post)', 0.2, 'r', time);

% Labels with larger font
xlabel('Time (ms)', 'FontSize', 14);
ylabel('Voltage (μV)', 'FontSize', 14);

% Set axis properties
ax = gca;
ax.YLim = [-30 30];
ax.YTick = -30:10:30;
ax.XLim = [0 550];
ax.XTick = 0:100:550;
ax.FontSize = 12;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'left';
ax.Box = 'off';
ax.TickDir = 'out';

% Tone onset windows [start, end] in ms
toneWindows = [0 116; 149 265; 298 414];

% Bar visual settings
barHeight = 4;
baseY = ax.YLim(1);        % bottom of y-axis
raisedY = baseY + 6;       % elevated bar height for 1st and 3rd tones

% Draw black bars
hold on;
for i = 1:size(toneWindows, 1)
    xBar = toneWindows(i, 1);
    width = toneWindows(i, 2) - toneWindows(i, 1);
    % First and third bars raised, second at bottom
    if i == 2
        yBar = baseY;
    else
        yBar = raisedY;
    end
    rectangle('Position', [xBar, yBar, width, barHeight], ...
              'FaceColor', 'k', 'EdgeColor', 'none');
end
hold off;

% Save figure (keep original filename)
saveas(gcf, '/home/ethan/fig_report/baseline_corrected_test_mean_ERP_triplet_response_for_channel_171_TYH.png');

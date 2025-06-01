% loads eeglab 
run('/home/ethan/scripts/eeglab2021.0/eeglab.m')

% load the data
load('/home/ethan/Data/BOD_volitional_streaming_04.mat')
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

%% plotting acouple of the channels (columns) as a sanity check without error bars

% Choose which column to plot 
columnIndex = 175;

% Extract the selected column
columnData = Triplet_baselinecorrected_Matrix_avrg(:, columnIndex);

% Create the x-axis (assuming time or row index)
x = 1:length(columnData);
time = x/srate;
% Plot the selected column
figure; % Create a new figure window
plot(time, columnData, 'b-', 'LineWidth', 2); 

% Add labels and title
xlabel('time (ms)'); % X-axis label
ylabel('voltage '); % Y-axis label
title(['Plot of Channel ', num2str(columnIndex)]); 

%% plot of all channels for comparison without error bars
% Get the number of channels
numChannels = size(Triplet_baselinecorrected_Matrix_avrg, 2);

% Define grid layout
gridRows = ceil(sqrt(numChannels)); 
gridCols = ceil(numChannels / gridRows);

% Define x-axis (assuming time or samples)
x = 1:size(Triplet_baselinecorrected_Matrix_avrg, 1);
time = x/srate*1000;
% Create figure
figure;

% Loop through each channel and create subplots in a grid layout
for i = 1:numChannels
    subplot(gridRows, gridCols, i); % Arrange in a grid
    plot(time, Triplet_baselinecorrected_Matrix_avrg(:, i), 'b', 'LineWidth', 1); % Plot channel in blue
    
    % Add x and y axes
    ax = gca;
    ax.XColor = 'k'; % Black x-axis
    ax.YColor = 'k'; % Black y-axis
    ax.XAxisLocation = 'origin'; % Force x-axis to go through y=0
    ax.YAxisLocation = 'left'; % Keep y-axis on the left side
    box off; % Remove surrounding box
    
    % Remove x and y ticks for cleaner look
    set(gca, 'XTick', [], 'YTick', []);

    % Ensure consistent y-limits across all plots
    yLimit = max(abs(Triplet_baselinecorrected_Matrix_avrg(:))); % Get max absolute value for symmetry
    ylim([-yLimit, yLimit]); % Make y-axis symmetric around 0
    
    % Add smaller channel number label in the bottom-right corner
    text(max(time) * 1, max(Triplet_baselinecorrected_Matrix_avrg(:)) * 0.9, num2str(i), ...
        'FontSize', 6, 'FontWeight', 'Bold', ...
        'Units', 'data', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');
    
    
end

% Add super title
sgtitle('Triplet Baseline corrected');
%saveas(gcf, '/home/ethan/fig/baseline_corrected_test_for_all_channels_without_errorbars.png');

%% Plot of single channel with error bars

% Choose which column to plot 
columnIndex = 175;

% Extract the selected column
columnData = Triplet_baselinecorrected_Matrix_avrg(:, columnIndex);

% Channels in first dimension
dp_data_ep = permute(Triplet_baselinecorrected_Matrix, [2 1 3]);

% Create the x-axis 
x = -60:length(columnData)-61;
time = x / srate * 1000;

% Compute symmetric y-axis limits
yMax = max(abs(columnData) + 5); % Find the largest absolute value
yLimits = [-yMax, yMax]; % Set symmetric limits

% Compute x-axis limits
xLimits = [min(time), max(time)]; % Ensure no unnecessary gaps

% Plot the selected column
figure; % Create a new figure window
stdshade_acj(squeeze(dp_data_ep(columnIndex,:,:))', 0.2, 'r', time);

% Add labels and title
xlabel('time (ms)'); % X-axis label
ylabel('voltage'); % Y-axis label
title(['Baseline Corrected Plot of Channel ', num2str(columnIndex)]); 

% Set axis properties
ax = gca; % Get current axes handle
ax.YLim = yLimits; % Apply symmetric y-axis limits
ax.XLim = xLimits; % Set x-axis limits to remove gaps
ax.XAxisLocation = 'origin'; % Move x-axis to y=0
ax.YAxisLocation = 'left'; % Keep y-axis on the left
ax.Box = 'off'; % Remove the surrounding box
ax.TickDir = 'out'; % Flip tick marks to the opposite side

% Add a vertical dashed line at x = 0 with transparency
hold on;
h = plot([0 0], yLimits, '--k', 'LineWidth', 1); % Dotted black line
h.Color = [0, 0, 0, 0.3]; % RGBA: Black (0,0,0) with 30% opacity
hold off;

% Save the figure
saveas(gcf, '/home/ethan/fig_BOD(p2)/baseline_corrected_test_for_single_channel_BOD(p2).png');

%% plot with SEM (error bars, standard error of the mean)

% Channels in first dimension
dp_data_ep = permute(Triplet_baselinecorrected_Matrix, [2 1 3]);

figure;

for chan = 1:size(dp_data_ep,1)
subplot(10,20,chan)
 stdshade_acj(squeeze(dp_data_ep(chan,:,:))',0.2,'r')

end

saveas(gcf, '/home/ethan/fig_BOD(p2)/baseline_corrected_test_for_all_channels_BOD(p2).png');

%% plot with SEM (error bars, standard error of the mean)

% channels in first dimension
dp_data_ep = permute(Triplet_baselinecorrected_Matrix, [2 1 3]);

numChannels = size(Triplet_baselinecorrected_Matrix_avrg, 2);

% Define grid layout
gridRows = ceil(sqrt(numChannels)); 
gridCols = ceil(numChannels / gridRows);

% Define x-axis (assuming time or samples)
x = 1:size(Triplet_baselinecorrected_Matrix_avrg, 1);
time = x/srate*1000;
% Create figure
figure;

% Loop through each channel and create subplots in a grid layout
for i = 1:numChannels
    subplot(gridRows, gridCols, i); % Arrange in a grid
        stdshade_acj(squeeze(dp_data_ep(i,:,:))',0.2,'r')
    
    % Add x and y axes
    ax = gca;
    ax.XColor = 'k'; % Black x-axis
    ax.YColor = 'k'; % Black y-axis
    ax.XAxisLocation = 'origin'; % Force x-axis to go through y=0
    ax.YAxisLocation = 'left'; % Keep y-axis on the left side
    box off; % Remove surrounding box
    
    % Remove x and y ticks for cleaner look
    set(gca, 'XTick', [], 'YTick', []);

    % Ensure consistent y-limits across all plots
    yLimit = max(abs(Triplet_baselinecorrected_Matrix_avrg(:))); % Get max absolute value for symmetry
    ylim([-yLimit, yLimit]); % Make y-axis symmetric around 0
    
    % Add smaller channel number label in the bottom-right corner
    text(max(time) * 1, max(Triplet_baselinecorrected_Matrix_avrg(:)) * 0.9, num2str(i), ...
        'FontSize', 6, 'FontWeight', 'Bold', ...
        'Units', 'data', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');
    
    
end

% Add super title
sgtitle('Triplet Baseline corrected');
saveas(gcf, '/home/ethan/fig_BOD(p2)/baseline_corrected_test_for_all_channels_BOD(p2).png');

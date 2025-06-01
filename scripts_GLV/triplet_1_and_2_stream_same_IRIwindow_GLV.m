clear
% loads eeglab 
run('/home/ethan/scripts/eeglab2021.0/eeglab.m')

% load the data 
load('/home/ethan/Data/GLV_volitional_streaming_04.mat')
whos

%% Finding the IRI (block size of 1s and 2s) for 1-stream and 2-stream perception

% Find indices where perception is 1
oneIndices = find(perception == 1);

% Identify breaks between separate sections of 1s
diffs1 = diff(oneIndices);
breaks1 = find(diffs1 > 1);

% Ensure startIndices and endIndices are defined properly for 1s
if isempty(breaks1)  % Case where all 1s are in one continuous block
    startIndices1 = oneIndices(1);
    endIndices1 = oneIndices(end);
else
    startIndices1 = [oneIndices(1); oneIndices(breaks1 + 1)];  
    endIndices1 = [oneIndices(breaks1); oneIndices(end)];
end

% Compute lengths of each 1s section
sectionLengths1 = endIndices1 - startIndices1 + 1;

% Find indices where perception is 2
twoIndices = find(perception == 2);

% Identify breaks between separate sections of 2s
diffs2 = diff(twoIndices);
breaks2 = find(diffs2 > 1);

% Ensure startIndices and endIndices are defined properly for 2s
if isempty(breaks2)  % Case where all 2s are in one continuous block
    startIndices2 = twoIndices(1);
    endIndices2 = twoIndices(end);
else
    startIndices2 = [twoIndices(1); twoIndices(breaks2 + 1)];  
    endIndices2 = [twoIndices(breaks2); twoIndices(end)];
end

% Compute lengths of each 2s section
sectionLengths2 = endIndices2 - startIndices2 + 1;

N1 = length(sectionLengths1);
N2 = length(sectionLengths2);
% Display results
disp('1 section lengths')
disp(sectionLengths1);
disp(['Number of 1 sections: ', num2str(length(sectionLengths1))]);
disp('2 section lengths')
disp(sectionLengths2);
disp(['Number of 2 sections: ', num2str(length(sectionLengths2))]);

% converting the lengths into time (S) for sampling rate 1200 Hz
IRI_1stream = sectionLengths1/1200;
IRI_2stream = sectionLengths2/1200;
disp(IRI_1stream);
disp(IRI_2stream);

%% calculating mean, median and mode of the IRI_1stream and IRI_2stream

% Compute statistics for blocks of 1s
meanOne = mean(IRI_1stream);
medianOne = median(IRI_1stream);
modeOne = mode(IRI_1stream);

% Compute statistics for blocks of 2s
meanTwo = mean(IRI_2stream);
medianTwo = median(IRI_2stream);
modeTwo = mode(IRI_2stream);

% Display results
fprintf('Statistics for blocks of 1s:\n');
fprintf('Mean: %.2f, Median: %.2f, Mode: %.2f\n', meanOne, medianOne, modeOne);

fprintf('\nStatistics for blocks of 2s:\n');
fprintf('Mean: %.2f, Median: %.2f, Mode: %.2f\n', meanTwo, medianTwo, modeTwo);

%% combining the two vectors to find the total IRI for both 1 and 2 streams 

% Combine the two vectors
combinedLengths = [IRI_1stream; IRI_2stream];

% Display statistics for the combined vector
meanCombined = mean(combinedLengths);
medianCombined = median(combinedLengths);
modeCombined = mode(combinedLengths);
numCombined = length(combinedLengths);

% Display the results
disp(['Mean: ', num2str(meanCombined)]);
disp(['Median: ', num2str(medianCombined)]);
disp(['Mode: ', num2str(modeCombined)]);
disp(['Total Count: ', num2str(numCombined)]);

%% converting the IRI median value (medianCombined) from seconds into "steps"  and round to an integer value 
medianIRI = round(medianCombined*1200);

disp(['medianIRI:', num2str(medianIRI)]);

%% Defining the windows of size medianIRI, centred around each value in 'startIndices1' and 'startIndices2'

% Compute window size
halfwindow = floor(medianIRI / 2);

% ensure that halfwindow is an integer value
halfWindow = round(halfwindow);

% Create windows around startIndices1 values
windowsOnes = arrayfun(@(idx) max(1, idx - halfWindow):min(length(perception), idx + halfWindow), startIndices1, 'UniformOutput', false);

% Create windows around startIndices2 values
windowsTwos = arrayfun(@(idx) max(1, idx - halfWindow):min(length(perception), idx + halfWindow), startIndices2, 'UniformOutput', false);

%% find the corresponding values in stim for these indices in the windows, ensuring output is a matrix and selected values retain their original indices that they had in stim

% Initialize output matrices with 0s (number of rows = number of windows, number of columns = length of stim)
stimOnesMatrix = zeros(length(windowsOnes), numel(stim));
stimTwosMatrix = zeros(length(windowsTwos), numel(stim));

% Assign extracted values while maintaining indices for windowsOnes
for i = 1:length(windowsOnes)
    stimOnesMatrix(i, windowsOnes{i}) = stim(windowsOnes{i});
end

% Assign extracted values while maintaining indices for windowsTwos
for i = 1:length(windowsTwos)
    stimTwosMatrix(i, windowsTwos{i}) = stim(windowsTwos{i});
end


%% identifying the indices at which the start of each tone occurs (i.e. the start of each block of 1s)

% Function to find first indices of 1s in separate blocks
findFirstOnesIndices = @(row) find(diff([0, row]) == 1);

% Apply function to each row of stimOnesMatrix
firstOnesIndicesOnes = cellfun(findFirstOnesIndices, num2cell(stimOnesMatrix, 2), 'UniformOutput', false);

% Apply function to each row of stimTwosMatrix
firstOnesIndicesTwos = cellfun(findFirstOnesIndices, num2cell(stimTwosMatrix, 2), 'UniformOutput', false);

%% Identifying the indices of the start of each triplet
% i.e. since each triplet contains 3 tones need to take every 3rd value
 
% subtracting 2 from the ends of each array otherwise it will include/count only part a triplet so subtracting 2 ensures that only full triplets are counted
% Apply to each row of firstOnesIndicesOnes
firstOnesIndicesOnesTrimmed = cellfun(@(x) x(1:min(end-2, end)), firstOnesIndicesOnes, 'UniformOutput', false);

% Apply to each row of firstOnesIndicesTwos 
firstOnesIndicesTwosTrimmed = cellfun(@(x) x(1:min(end-2, end)), firstOnesIndicesTwos, 'UniformOutput', false);

% Function to extract every 3rd value starting from the first value
extractEveryThird = @(arr) arr(1:3:end);

% Apply function to each row of firstOnesIndicesOnesTrimmed
everyThirdOnes = cellfun(extractEveryThird, firstOnesIndicesOnesTrimmed, 'UniformOutput', false);

% Apply function to each row of firstOnesIndicesTwosTrimmed
everyThirdTwos = cellfun(extractEveryThird, firstOnesIndicesTwosTrimmed, 'UniformOutput', false);

% Convert cell arrays into separate vectors
vectorEveryThirdOnes = [everyThirdOnes{:}]; % Flatten into a single vector
vectorEveryThirdTwos = [everyThirdTwos{:}]; % Flatten into a single vector

% Rename the vectors
selectedIndices_stimOne = vectorEveryThirdOnes;  
selectedIndices_stimTwo = vectorEveryThirdTwos;  

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

%% Plot of single channel with error bars

% Choose which column to plot 
columnIndex = 2;

% Extract the selected column from both matrices
columnData_1stream = Triplet_baselinecorrected_Matrix_1stream(:, columnIndex);
columnData_2streams = Triplet_baselinecorrected_Matrix_2streams(:, columnIndex);

% Channels in first dimension for both matrices
dp_data_ep_1stream = permute(Triplet_baselinecorrected_Matrix_1stream, [2 1 3]);
dp_data_ep_2streams = permute(Triplet_baselinecorrected_Matrix_2streams, [2 1 3]);

% Create the x-axis 
x = -60:length(columnData_1stream)-61;
time = x / srate * 1000;

% Compute symmetric y-axis limits
yMax = max([max(abs(columnData_1stream)), max(abs(columnData_2streams))]) + 25; % Find the largest absolute value
yLimits = [-500, 500]; % Set symmetric limits

% Compute x-axis limits
xLimits = [min(time), max(time)]; % Ensure no unnecessary gaps

% Plot the selected column for both datasets
figure; % Create a new figure window

hold on;
stdshade_acj(squeeze(dp_data_ep_1stream(columnIndex,:,:))', 0.2, 'r', time); % Plot 1-stream in red
stdshade_acj(squeeze(dp_data_ep_2streams(columnIndex,:,:))', 0.2, 'b', time); % Plot 2-streams in blue
hold off;

% Add labels and title
xlabel('Time (ms)'); % X-axis label
ylabel('Voltage'); % Y-axis label
title(['Baseline and IRI Median Windowed Plot of Channel ', num2str(columnIndex)]); 

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

% Manually add legend with invisible lines
h1 = plot(NaN, NaN, 'r', 'LineWidth', 2); % Red line for 1-stream
h2 = plot(NaN, NaN, 'b', 'LineWidth', 2); % Blue line for 2-streams
legend([h1, h2], {'1-stream', '2-streams'}, 'Location', 'northeast');

hold off;

% Save the figure
saveas(gcf, '/home/ethan/fig_GLV(p3)/baseline_corrected_and_IRIMedian_combined_windowed_separated_triplets_for_channel2_GLV(p3).png');

%% Plot with SEM (Standard Error of the Mean)

% Channels in first dimension for both matrices
dp_data_ep_1stream = permute(Triplet_baselinecorrected_Matrix_1stream, [2 1 3]);
dp_data_ep_2streams = permute(Triplet_baselinecorrected_Matrix_2streams, [2 1 3]);

numChannels = size(Triplet_baselinecorrected_Matrix_1stream, 2);

% Define grid layout
gridRows = ceil(sqrt(numChannels)); 
gridCols = ceil(numChannels / gridRows);

% Define x-axis (assuming time or samples)
x = 1:size(Triplet_baselinecorrected_Matrix_1stream, 1);
time = x / srate * 1000;

% Create figure
figure;

% Loop through each channel and create subplots in a grid layout
for i = 1:numChannels
    subplot(gridRows, gridCols, i); % Arrange in a grid
    hold on;
    
    % Plot both 1-stream (red) and 2-streams (blue)
    stdshade_acj(squeeze(dp_data_ep_1stream(i,:,:))', 0.2, 'r', time);
    stdshade_acj(squeeze(dp_data_ep_2streams(i,:,:))', 0.2, 'b', time);
    
    hold off;
    
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
    yLimit = max(abs([Triplet_baselinecorrected_Matrix_1stream(:); Triplet_baselinecorrected_Matrix_2streams(:)])+70); % Get max absolute value for symmetry
    ylim([-25,25]); % Make y-axis symmetric around 0
    
    % Add smaller channel number label in the bottom-right corner
    text(max(time) * 1, 25 * 0.9, num2str(i), ...
        'FontSize', 6, 'FontWeight', 'Bold', ...
        'Units', 'data', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');
    
end

% Add super title
sgtitle('1 stream and 2 stream Baseline and IRI Median Windowed Triplets for All Channels');

% Save the figure
saveas(gcf, '/home/ethan/fig_GLV(p3)/baseline_corrected_and_IRIMedian_combined_windowed_separated_triplets_for_all_channels_GLV(p3).png');



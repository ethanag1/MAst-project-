clear 
clc
%% setup
dir_coi = '/home/ethan/coi'; % set filepath to coi toolbox on pc
cd(dir_coi) % change matlab's current directory to it
coi_init % make coi code available for matlab to use
% make sure to also add the gcmi toolbox to path
addpath('/home/ethan/gcmi/matlab');
savepath

% loads eeglab 
run('/home/ethan/scripts/eeglab2021.0/eeglab.m')

% load the data 
load('/home/ethan/Data/BOD_volitional_streaming_04.mat')

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

%% combining the two types of triplets identified into 1 single vector retaining the correct trial ordering
% Ensure both vectors are column vectors
selectedIndices_stimOne = selectedIndices_stimOne(:);  
selectedIndices_stimTwo = selectedIndices_stimTwo(:);  

% Now combine them
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
BaselineMatrix_3D_trimmed = Triplet_baselinecorrected_Matrix(61:end, :, :);

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
    n_t = 661; % number of timepoints
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
t_ = 0:1:660;
n_t = min(n_t, numel(t_)); % this stops selection of more timepoints than there are in the input
t_ = t_(1:n_t)/1.2; % convert to seconds

%% calculating the mi and zmi values for each channel 
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

% Display results
disp('Significant channels:');
disp(significant_channels);
disp('number of channels with significant MI values');
disp(length(significant_channels));

%% ii calculation for significant channels 

n_pl3 = 100; % number of permutations for COI
% Initialize storage for ii and zii values (as cell arrays)
ii_all = cell(1, num_channels);
zii_all = cell(1, num_channels);


% Loop over each significant channel
for ch = significant_channels
    % Extract data for the current significant channel
    rr_ = BaselineMatrix_3D_ep(:, :, ch);
    
    % Compute ii and its statistical significance
    [ii_, zii_] = ii_decode(l_, rr_, 't', n_t,'pl', n_pl3);
    % Store results in cell arrays
    ii_all{ch} = ii_;  % Store full 2D matrix
    zii_all{ch} = zii_;
    
    
    % Display progress
    %fprintf('Processed Channel %d/%d\n', find(significant_channels == ch), length(significant_channels));
end

% Output confirmation
disp('✅ ii computed for all significant channels.');


%% Plotting the significant COI values for all significant channels

% Define grid layout based on the number of significant channels
num_sig_channels = length(significant_channels);
num_rows = ceil(sqrt(num_sig_channels));  % Arrange channels in a square grid
num_cols = ceil(num_sig_channels / num_rows);

% Create a figure
figure;
clf;

% Loop over each significant channel and plot its COI heatmap
for i = 1:num_sig_channels
    ch = significant_channels(i); % Get channel index

    % Check if it's the last subplot (to make it larger)
    if i == num_sig_channels
        % Use the position property to enlarge the last subplot (make room for colorbar)
        ax = subplot(num_rows, num_cols, i);
        pos = get(ax, 'Position');
        pos(3) = pos(3) * 1.43; % Increase width (or adjust as needed)
        pos(4) = pos(4) * 0.96; % Increase height (or adjust as needed)
        set(ax, 'Position', pos); % Apply the new position
    else
        % Regular subplot for all but the last
        subplot(num_rows, num_cols, i);
    end
    
    % Extract COI values (now stored in a cell array)
    zii_data = zii_all{ch};  % ✅ Corrected indexing

    % Extract the COI field from zii_data
    zcoi_data = zii_data.coi;

    % Plot COI heatmap for the current channel
    z_heatmap(zcoi_data, t_); 
    
    % Adjust color axis to the range [0, 1]
    clim([0, 1]);

    % Remove tick labels and axis labels
    
    % Remove tick labels and axis labels (for all subplots except the last one)
    set(gca, 'XTick', [], 'YTick', [], 'YTickLabel', [], 'FontSize', 1);  % Set Y-axis font size to 1

    box on; % Keeps the box around each subplot
    % Add channel number in the title
    title([ num2str(ch)], 'FontSize', 7, 'FontWeight', 'bold');

   
    if i == num_sig_channels
        cb = colorbar; % Add the colorbar
        % Set colorbar ticks to only display the min and max values
        set(cb, 'Limits', [0, 1]);  % Set colorbar limits directly
        cb.Ticks = [0, 1]; % Only display min and max ticks on the colorbar
        % Only for the last subplot, show the axis labels with a larger font size
        xlabel('Time (ms)', 'FontSize', 7);  % Larger X-axis label
        ylabel('Time (ms)', 'FontSize', 7);  % Larger Y-axis label
        % Add a label to the colorbar
        ylabel(cb, 'Signif.CoI {(z-score)}', 'FontSize', 7); % Colorbar label with font size
        set(cb, 'FontSize', 7);  % Increase font size of the ticks
    else
        colorbar off;  % Remove colorbar for other subplots
    end

   
end

% Add global title
sgtitle(' Significant Co-I Values for Significant Channels IRI Window method 100 Permutations');

% Save the figure
saveas(gcf, '/home/ethan/fig_BOD(p2)/within_channel_signifcoi_significant_channels_same_IRIwindow_BOD(p2).png');
disp('✅ Saved COI plot for all significant channels.');

%% Plotting the significant synergistic COI values for all significant channels

% Define grid layout based on the number of significant channels
num_sig_channels = length(significant_channels);
num_rows = ceil(sqrt(num_sig_channels));  % Arrange channels in a square grid
num_cols = ceil(num_sig_channels / num_rows);

% Create a figure
figure;
clf;

% Loop over each significant channel and plot its COI heatmap
for i = 1:num_sig_channels
    ch = significant_channels(i); % Get channel index

    % Check if it's the last subplot (to make it larger)
    if i == num_sig_channels
        % Use the position property to enlarge the last subplot (make room for colorbar)
        ax = subplot(num_rows, num_cols, i);
        pos = get(ax, 'Position');
        pos(3) = pos(3) * 1.43; % Increase width (or adjust as needed)
        pos(4) = pos(4) * 0.96; % Increase height (or adjust as needed)
        set(ax, 'Position', pos); % Apply the new position
    else
        % Regular subplot for all but the last
        subplot(num_rows, num_cols, i);
    end
    
    % Extract COI values (now stored in a cell array)
    zii_data = zii_all{ch};  % ✅ Corrected indexing


    % Extract the COI field from zii_data
    zcoi_data = abs(zii_data.syn);

    % Plot COI heatmap for the current channel
    z_heatmap(zcoi_data, t_); 
    
    % Adjust color axis to the range [0, 1]
    clim([0, 1]);

    % Remove tick labels and axis labels
    
    % Remove tick labels and axis labels (for all subplots except the last one)
    set(gca, 'XTick', [], 'YTick', [], 'YTickLabel', [], 'FontSize', 1);  % Set Y-axis font size to 1

    box on; % Keeps the box around each subplot
    % Add channel number in the title
    title([ num2str(ch)], 'FontSize', 7, 'FontWeight', 'bold');

   
    if i == num_sig_channels
        cb = colorbar; % Add the colorbar
        % Set colorbar ticks to only display the min and max values
        set(cb, 'Limits', [0, 1]);  % Set colorbar limits directly
        cb.Ticks = [0, 1]; % Only display min and max ticks on the colorbar
        % Only for the last subplot, show the axis labels with a larger font size
        xlabel('Time (ms)', 'FontSize', 7);  % Larger X-axis label
        ylabel('Time (ms)', 'FontSize', 7);  % Larger Y-axis label
        % Add a label to the colorbar
        ylabel(cb, 'Signif.CoI {(z-score)}', 'FontSize', 7); % Colorbar label with font size
        set(cb, 'FontSize', 7);  % Increase font size of the ticks
    else
        colorbar off;  % Remove colorbar for other subplots
    end

   
end

% Add global title
sgtitle(' Significant Synergistic Co-I Values for Significant Channels IRI Window method 100 Permutations', 'FontSize', 10);

% Save the figure
saveas(gcf, '/home/ethan/fig_BOD(p2)/within_channel_syn_signifcoi_significant_channels_same_IRIwindow_BOD(p2).png');
disp('✅ Saved COI plot for all significant channels.');

%% Plotting the significant redundant COI values for all significant channels

% Define grid layout based on the number of significant channels
num_sig_channels = length(significant_channels);
num_rows = ceil(sqrt(num_sig_channels));  % Arrange channels in a square grid
num_cols = ceil(num_sig_channels / num_rows);

% Create a figure
figure;
clf;

% Loop over each significant channel and plot its COI heatmap
for i = 1:num_sig_channels
    ch = significant_channels(i); % Get channel index

    % Check if it's the last subplot (to make it larger)
    if i == num_sig_channels
        % Use the position property to enlarge the last subplot (make room for colorbar)
        ax = subplot(num_rows, num_cols, i);
        pos = get(ax, 'Position');
        pos(3) = pos(3) * 1.43; % Increase width (or adjust as needed)
        pos(4) = pos(4) * 0.96; % Increase height (or adjust as needed)
        set(ax, 'Position', pos); % Apply the new position
    else
        % Regular subplot for all but the last
        subplot(num_rows, num_cols, i);
    end
    
    % Extract COI values (now stored in a cell array)
    zii_data = zii_all{ch};  % ✅ Corrected indexing

    % Extract the COI field from zii_data
    zcoi_data = zii_data.red;

    % Plot COI heatmap for the current channel
    z_heatmap(zcoi_data, t_); 
    
    % Adjust color axis to the range [0, 1]
    clim([0, 1]);

    % Remove tick labels and axis labels
    
    % Remove tick labels and axis labels (for all subplots except the last one)
    set(gca, 'XTick', [], 'YTick', [], 'YTickLabel', [], 'FontSize', 1);  % Set Y-axis font size to 1

    box on; % Keeps the box around each subplot
    % Add channel number in the title
    title([ num2str(ch)], 'FontSize', 7, 'FontWeight', 'bold');

   
    if i == num_sig_channels
        cb = colorbar; % Add the colorbar
        % Set colorbar ticks to only display the min and max values
        set(cb, 'Limits', [0, 1]);  % Set colorbar limits directly
        cb.Ticks = [0, 1]; % Only display min and max ticks on the colorbar
        % Only for the last subplot, show the axis labels with a larger font size
        xlabel('Time (ms)', 'FontSize', 7);  % Larger X-axis label
        ylabel('Time (ms)', 'FontSize', 7);  % Larger Y-axis label
        % Add a label to the colorbar
        ylabel(cb, 'Signif.CoI {(z-score)}', 'FontSize', 7); % Colorbar label with font size
        set(cb, 'FontSize', 7);  % Increase font size of the ticks
    else
        colorbar off;  % Remove colorbar for other subplots
    end

   
end

% Add global title
sgtitle('Significant Redundant Co-I Values for Significant Channels IRI Window method 100 Permutations', 'FontSize', 10);

% Save the figure
saveas(gcf, '/home/ethan/fig_BOD(p2)/within_channel_red_signifcoi_significant_channels_same_IRIwindow_BOD(p2).png');
disp('✅ Saved COI plot for all significant channels.');

%% Plotting the COI values for all significant channels

% Define grid layout based on the number of significant channels
num_sig_channels = length(significant_channels);
num_rows = ceil(sqrt(num_sig_channels));  % Arrange channels in a square grid
num_cols = ceil(num_sig_channels / num_rows);

% Create a figure
figure;
clf;

% Loop over each significant channel and plot its COI heatmap
for i = 1:num_sig_channels
    ch = significant_channels(i); % Get channel index

    % Check if it's the last subplot (to make it larger)
    if i == num_sig_channels
        % Use the position property to enlarge the last subplot (make room for colorbar)
        ax = subplot(num_rows, num_cols, i);
        pos = get(ax, 'Position');
        pos(3) = pos(3) * 1.43; % Increase width (or adjust as needed)
        pos(4) = pos(4) * 0.96; % Increase height (or adjust as needed)
        set(ax, 'Position', pos); % Apply the new position
    else
        % Regular subplot for all but the last
        subplot(num_rows, num_cols, i);
    end
    
    % Extract COI values (now stored in a cell array)
    coi_ = ii_all{ch};  % ✅ Corrected indexing
    coi_data = coi_.coi; % Extract the COI field from coi_data
    % Plot COI heatmap for the current channel
    coi_heatmap(coi_data, t_); 
    
    % Adjust color axis to the range [-0.1, 0.1]
    clim([-0.01, 0.01]);

    % Remove tick labels and axis labels
    
    % Remove tick labels and axis labels (for all subplots except the last one)
    set(gca, 'XTick', [], 'YTick', [], 'YTickLabel', [], 'FontSize', 1);  % Set Y-axis font size to 1

    box on; % Keeps the box around each subplot
    % Add channel number in the title
    title([ num2str(ch)], 'FontSize', 7, 'FontWeight', 'bold');

   
    if i == num_sig_channels
        cb = colorbar; % Add the colorbar
        % Set colorbar ticks to only display the min and max values
        clim([-0.01 0.01]);  % Set the color axis limits (same as clim)
        cb.Ticks = [-0.01, 0.01]; % Only display min and max ticks on the colorbar
        
        % Only for the last subplot, show the axis labels with a larger font size
        xlabel('Time (ms)', 'FontSize', 7);  % Larger X-axis label
        ylabel('Time (ms)', 'FontSize', 7);  % Larger Y-axis label
        % Add a label to the colorbar
        ylabel(cb, 'CoI_{log_{2}n}', 'FontSize', 7); % Colorbar label with font size
        set(cb, 'FontSize', 7);  % Increase font size of the ticks
    else
        colorbar off;  % Remove colorbar for other subplots
    end

   
end

% Add global title
sgtitle('Co-Information Heatmaps for Significant Channels IRI Window method');

% Save the figure
saveas(gcf, '/home/ethan/fig_BOD(p2)/within_channel_coi_significant_channels_same_IRIwindow_BOD(p2).png');
disp('✅ Saved COI plot for all significant channels.');

%% Display results
disp('Significant channels:');
disp(significant_channels);
disp('number of channels with significant MI values');
disp(length(significant_channels));

%% percentage of triplets identified using exclusion method
Percentage_triplets_identified = length(combinedIndices) / 500 * 100;
disp('Percentage of triplets retained using direct extraction method');
disp(Percentage_triplets_identified);
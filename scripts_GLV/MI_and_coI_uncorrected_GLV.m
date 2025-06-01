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
load('/home/ethan/Data/GLV_volitional_streaming_04.mat')

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

%%
% Define number of channels and timepoints
[num_timepoints, num_channels] = size(mi_all);

% Define grid layout
num_rows = ceil(sqrt(num_channels));  % Arrange channels in a square grid
num_cols = ceil(num_channels / num_rows);

% Determine global y-axis limits
y_min = min(mi_all(:));
y_max = max(mi_all(:)) * 1.2; % Set a slightly higher max for uniformity
sig_bar_y = max(mi_all(:)) * 1.05; % Fixed position for significance bars

figure; % Create a new figure
for ch = 1:num_channels
    subplot(num_rows, num_cols, ch); % Create subplot for each channel
    hold on;
    
    % Plot MI values in black
    plot(1:num_timepoints, mi_all(:, ch), 'k', 'LineWidth', 1);
    
    % Find significant time points where zmi_ > 1.645
    significant_points = find(zmi_all(:, ch) > 1.645);
    
    % Plot purple bars at the same y position across all plots
    if ~isempty(significant_points)
        for t = significant_points'
            plot([t t], [sig_bar_y, sig_bar_y * 1.1], 'Color', [0.5 0 0.5], 'LineWidth', 1);
        end
    end
    
    % Add channel number in the top-right in smaller font
    text(num_timepoints * 1, y_max * 1, num2str(ch), ...
        'FontSize', 5, 'FontWeight', 'bold', 'HorizontalAlignment', 'right');

    % Set consistent y-axis limits for all subplots
    ylim([y_min y_max]);

    % Set X-axis limits to remove any gap on the right
    xlim([1, num_timepoints]);

    % Remove subplot box
    set(gca, 'Box', 'off');

    % Remove x and y ticks for a cleaner look
    xticks([]); yticks([]);
    
    hold off;
end


sgtitle('Mutual Information for All Channels with Significant Time Points'); % Add a global title

saveas(gcf, '/home/ethan/fig_GLV(p3)/mi_all_channels_significant_uncorrected_GLV(p3).png');

 %% Log plotting the MI values

% Select the channel to plot
channel_to_plot = 1; % Change this to plot a different channel

% Define timepoints in milliseconds
num_timepoints = size(mi_all, 1);
sampling_rate = 1200; % Hz
time_ms = (1:num_timepoints) / sampling_rate * 1000; % Convert time to milliseconds

figure; % Create new figure
hold on;

% Plot MI values in black with log base 2 scale
semilogy(time_ms, mi_all(:, channel_to_plot), 'k', 'LineWidth', 1.5);
set(gca, 'YScale', 'log'); % Set logarithmic scale

% Adjust log base to 2 by modifying YTick labels
yticks(2.^floor(log2(min(mi_all(:, channel_to_plot)))) : 2.^ceil(log2(max(mi_all(:, channel_to_plot)))));
ylabels = arrayfun(@(x) sprintf('2^{%d}', log2(x)), yticks, 'UniformOutput', false);
set(gca, 'YTickLabel', ylabels);

% Find significant time points where zmi_ > 1.645
significant_points = find(zmi_all(:, channel_to_plot) > 1.645);

% Plot purple bars for significant time points
if ~isempty(significant_points)
    y_max = max(mi_all(:, channel_to_plot)) * 1.1; % Slightly above max MI
    for t = significant_points'
        plot([time_ms(t) time_ms(t)], [y_max, y_max * 1.15], 'Color', [0.5 0 0.5], 'LineWidth', 1.5);
    end
end

% Add labels and title
xlabel('Time (ms)');
ylabel('MI_{log_{2}}');
title(['Mutual Information for Channel ', num2str(channel_to_plot)]);

% Customize x-axis ticks to increment by 50 ms
xticks(0:50:max(time_ms));

% Set X-axis limits to remove any gap on the right
xlim([min(time_ms), max(time_ms)]);

% Remove box around the plot
set(gca, 'Box', 'off');

hold off;

% Save figure
saveas(gcf, '/home/ethan/fig_GLV(p3)/mi_channel1_significantlog_uncorrected_GLV(p3).png');


%% normal scale plotting the mi values
% Select the channel to plot
channel_to_plot = 1; % Change this to plot a different channel

% Define timepoints in milliseconds
num_timepoints = size(mi_all, 1);
sampling_rate = 1200; % Hz
time_ms = (1:num_timepoints) / sampling_rate * 1000; % Convert to ms

figure; % Create new figure
hold on;

% Plot MI values in black (normal Y-axis)
plot(time_ms, mi_all(:, channel_to_plot), 'k', 'LineWidth', 1.5);

% Find significant time points where zmi_ > 1.645
significant_points = find(zmi_all(:, channel_to_plot) > 1.645);

% Plot purple bars for significant time points
if ~isempty(significant_points)
    y_max = max(mi_all(:, channel_to_plot)) * 1.1; % Slightly above max MI
    for t = significant_points'
        plot([time_ms(t) time_ms(t)], [y_max, y_max * 1.05], 'Color', [0.5 0 0.5], 'LineWidth', 1.1);
    end
end

% Add labels and title
xlabel('Time (ms)');
ylabel('Mutual Information');
title(['Mutual Information for Channel ', num2str(channel_to_plot)]);

% Customize x-axis ticks to increment by 50 ms
xticks(0:50:max(time_ms));

% Set X-axis limits to remove any gap on the right
xlim([min(time_ms), max(time_ms)]);

% Remove box around the plot
set(gca, 'Box', 'off');

hold off;
saveas(gcf, '/home/ethan/fig_GLV(p3)/mi_channel1_significantnormal_uncorrected_GLV(p3).png');
disp(max(time_ms));

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

%% calculating the within-channel jmi and coi values for each channel with significant MI values

n_pl2 = 100; % number of permutations for COI
% Initialize storage for JMI and COI values (as cell arrays)
jmi_all = cell(1, num_channels);
zjmi_all = cell(1, num_channels);
coi_all = cell(1, num_channels);

% Loop over each significant channel
for ch = significant_channels
    % Extract data for the current significant channel
    rr_ = BaselineMatrix_3D_ep(:, :, ch);
    
    % Compute Joint Mutual Information (JMI) and its statistical significance
    [jmi_, zjmi_] = jmi_decode(l_, rr_, 't', n_t, 'pl', n_pl2);
    
    % Compute Co-Information (COI) as the difference between MI and JMI
    coi_ = coi(mi_all(:, ch), jmi_);
    
    % Store results in cell arrays
    jmi_all{ch} = jmi_;  % Store full 2D matrix
    zjmi_all{ch} = zjmi_;
    coi_all{ch} = coi_;
    
    
    % Display progress
    fprintf('Processed Channel %d/%d\n', find(significant_channels == ch), length(significant_channels));
end

% Output confirmation
disp('✅ ii computed for all significant channels.');

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
        pos(3) = pos(3) * 1.5; % Increase width (or adjust as needed)
        pos(4) = pos(4) * 1; % Increase height (or adjust as needed)
        set(ax, 'Position', pos); % Apply the new position
    else
        % Regular subplot for all but the last
        subplot(num_rows, num_cols, i);
    end
    
    % Extract COI values (now stored in a cell array)
    coi_data = coi_all{ch};  % ✅ Corrected indexing

    % Plot COI heatmap for the current channel
    coi_heatmap(coi_data, t_); 
    
    % Adjust color axis to the range [-0.1, 0.1]
    clim([-0.01, 0.01]);

    % Remove tick labels and axis labels
    
    % Remove tick labels and axis labels (for all subplots except the last one)
    set(gca, 'XTick', [], 'YTick', [], 'YTickLabel', [], 'FontSize', 1);  % Set Y-axis font size to 1

    box on; % Keeps the box around each subplot
    % Add channel number in the title
    title([ num2str(ch)], 'FontSize', 5, 'FontWeight', 'bold');

   
    if i == num_sig_channels
        cb = colorbar; % Add the colorbar
        % Set colorbar ticks to only display the min and max values
        clim([-0.01 0.01]);  % Set the color axis limits (same as clim)
        cb.Ticks = [-0.01, 0.01]; % Only display min and max ticks on the colorbar
        
        % Only for the last subplot, show the axis labels with a larger font size
        xlabel('Time (ms)', 'FontSize', 6);  % Larger X-axis label
        ylabel('Time (ms)', 'FontSize', 6);  % Larger Y-axis label
        % Add a label to the colorbar
        ylabel(cb, 'CoI_{log_{2}n}', 'FontSize', 6); % Colorbar label with font size
        set(cb, 'FontSize', 6);  % Increase font size of the ticks
    else
        colorbar off;  % Remove colorbar for other subplots
    end

   
end

% Add global title
sgtitle('Co-Information Heatmaps for Significant Channels Uncorrected method');

% Save the figure
saveas(gcf, '/home/ethan/fig_GLV(p3)/within_channel_coi_significant_channels_uncorrected_GLV(p3).png');
disp('✅ Saved COI plot for all significant channels.');

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
    fprintf('Processed Channel %d/%d\n', find(significant_channels == ch), length(significant_channels));
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
        pos(3) = pos(3) * 1.5; % Increase width (or adjust as needed)
        pos(4) = pos(4) * 1; % Increase height (or adjust as needed)
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
    title([ num2str(ch)], 'FontSize', 5, 'FontWeight', 'bold');

   
    if i == num_sig_channels
        cb = colorbar; % Add the colorbar
        % Set colorbar ticks to only display the min and max values
        set(cb, 'Limits', [0, 1]);  % Set colorbar limits directly
        cb.Ticks = [0, 1]; % Only display min and max ticks on the colorbar
        % Only for the last subplot, show the axis labels with a larger font size
        xlabel('Time (ms)', 'FontSize', 6);  % Larger X-axis label
        ylabel('Time (ms)', 'FontSize', 6);  % Larger Y-axis label
        % Add a label to the colorbar
        ylabel(cb, 'Signif.CoI {(z-score)}', 'FontSize', 6); % Colorbar label with font size
        set(cb, 'FontSize', 6);  % Increase font size of the ticks
    else
        colorbar off;  % Remove colorbar for other subplots
    end

   
end

% Add global title
sgtitle(' Significant Co-I Values for Significant Channels Uncorrected method');

% Save the figure
saveas(gcf, '/home/ethan/fig_GLV(p3)/within_channel_signifcoi_significant_channels_uncorrected_GLV(p3).png');
disp('✅ Saved COI plot for all significant channels.');

%% Generate and save separate COI plots for each significant channel

% Define output directory for saving figures
output_dir = '/home/ethan/fig_GLV(p3)/co-I for significant channels uncorrected method';  
if ~exist(output_dir, 'dir')
    mkdir(output_dir); % Create the directory if it doesn't exist
end

% Loop over each significant channel
for i = 1:length(significant_channels)
    ch = significant_channels(i); % Get channel index
    
    % Extract COI values for the current channel
    coi_data = coi_all{ch};  % ✅ Corrected indexing

    % Create a new figure
    figure;
    clf;
    
    % Plot COI heatmap for the current channel
    coi_heatmap(coi_data, t_); 

    % Adjust color axis to the range [-0.1, 0.1]
    clim([-0.01, 0.01]);

    % Modify colorbar and ensure it matches the new limits
    h = colorbar; 
    set(h, 'Limits', [-0.01, 0.01]);  % Explicitly set colorbar range
    ylabel(h, 'CoI_{log_{2}n}'); % Add label to colorbar

    % Add title with channel number
    title(['Channel ', num2str(ch)], 'FontSize', 10, 'FontWeight', 'bold');

    % Define filename for saving the figure
    filename = fullfile(output_dir, ['within_channel_coi_channel_' num2str(ch) '_uncorrected_GLV(p3).png']);
    
    % Save the figure
    saveas(gcf, filename);

    % Close figure to free up memory
    close(gcf);  
    
    % Display progress
    fprintf('✅ Saved COI plot for Channel %d: %s\n', ch, filename);
end

%% percentage of triplets identified using exclusion method
Percentage_triplets_identified = length(combinedIndices) / 500 * 100;
disp('Percentage of triplets retained using direct extraction method');
disp(Percentage_triplets_identified);

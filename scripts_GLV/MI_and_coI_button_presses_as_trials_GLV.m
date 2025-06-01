clear 
clc
%% Initialize parallel pool with up to 24 workers
%if isempty(gcp('nocreate'))
 %   parpool('local', 24);  % Adjust based on server's available cores
%end
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

%% Find start indices of each block of 1s and 2s separately

% Find where perception is 1
oneIndices = find(perception == 1);
diffs1 = diff(oneIndices);
breaks1 = find(diffs1 > 1);

% First index of each block of 1s
if isempty(breaks1)
    firstOnes = oneIndices(1);
else
    firstOnes = [oneIndices(1); oneIndices(breaks1 + 1)];
end

% Find where perception is 2
twoIndices = find(perception == 2);
diffs2 = diff(twoIndices);
breaks2 = find(diffs2 > 1);

% First index of each block of 2s
if isempty(breaks2)
    firstTwos = twoIndices(1);
else
    firstTwos = [twoIndices(1); twoIndices(breaks2 + 1)];
end

% renaming the variables
selectedIndices_stimOne = firstOnes; % Indices for 2 stream to 1 stream button presses
selectedIndices_stimTwo = firstTwos; % Indices for 1 stream to 2 stream button presses

%% combining the two types of triplets identified into 1 single vector retaining the correct trial ordering
% Ensure both vectors are column vectors
selectedIndices_stimOne = selectedIndices_stimOne(:); % Indices for 2 stream to 1 stream button presses
selectedIndices_stimTwo = selectedIndices_stimTwo(:); % Indices for 1 stream to 2 stream button presses 

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

%% creating sections of length -1200 and +1200 entries around each selected index 
% (550ms i.e. the length of the triplet) for each entry in "selectedIndices" and making cuts 
% in dp_data and transforming the 2D dp_data matrix into a 3d one with the
% 3rd dim being each section cut from the dp_data matrix

% Define the range around each selected index
preTime = 1200;    % Number of points before the selected index (baseline period)
postTime = 1200;  % Number of points after the selected index (triplet period)
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



%% permuting the inputs so that the matrix is trials x samples x channels

BaselineMatrix_3D_ep = permute(BaselineMatrix_3D, [3 1 2]);

% Display confirmation
disp('✅ Successfully permuted matrix.');
disp(['Size of 3D matrix: ', num2str(size(BaselineMatrix_3D_ep,1)), ' x ', num2str(size(BaselineMatrix_3D_ep,2)), ' x ', num2str(size(BaselineMatrix_3D_ep,3))]);

%% testing parameters
% edit these to change:
    % how many timepoints you compute mi/coi on
    % how many permutations you run for significance testing
    n_t = 2401; % number of timepoints
    n_pl = 100; % number of permutations
%% creating vector that labels the trials in the input according to the stimulus

% Initialize l_ vector with zeros
l_ = zeros(size(combinedIndices));

% Label trials: 0 for two to one stream, 1 for one to two stream
l_(ismember(combinedIndices, selectedIndices_stimOne)) = 0;
l_(ismember(combinedIndices, selectedIndices_stimTwo)) = 1;
disp('length of l_ vector created');
disp(length(l_));
disp(l_);
%% extract inputs

% to neatly plot the mi, we also need the vector labelling timepoints
t_ = -1200:1:1200;
n_t = min(n_t, numel(t_)); % this stops selection of more timepoints than there are in the input
t_ = (-1200:1200) / 1.2; % convert to ms

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
    [ii_, zii_] = ii_decode_edited(l_, rr_, 't', n_t, 'pl', n_pl3, 'is_x', false);

    
    % Store results directly into the main cell arrays
    ii_all{ch} = ii_;
    zii_all{ch} = zii_;
end



% Loop over each significant channel
%parfor idx = 1:length(significant_channels)
    %ch = significant_channels(idx);
    % Extract data for the current significant channel
    %rr_ = BaselineMatrix_3D_ep(:, :, ch);
    
    % Compute ii and its statistical significance
    %[ii_, zii_] = ii_decode_edited(l_, rr_, 't', n_t,'pl', n_pl3);
    
    % Use temporary variables to store results
    %temp_ii{idx} = ii_;  % Store full 2D matrix
    %temp_zii{idx} = zii_;
%end

% Assign results back to the main cell arrays
%for idx = 1:length(significant_channels)
    %ch = significant_channels(idx);
    %ii_all{ch} = temp_ii{idx};
    %zii_all{ch} = temp_zii{idx};
%end

% Output confirmation
disp('✅ ii computed for all significant channels.');

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

saveas(gcf, '/home/ethan/fig_GLV(p3)/mi_all_channels_significant_buttonpress_trials_GLV(p3).png');

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
sgtitle(' Significant Co-I Values for Significant Channels Button presses as trials 100 Permutations');

% Save the figure
saveas(gcf, '/home/ethan/fig_GLV(p3)/within_channel_signifcoi_significant_channels_buttonpress_trials_GLV(p3).png');
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
sgtitle(' Significant Synergistic Co-I Values for Significant Channels Button presses as trials 100 Permutations', 'FontSize', 10);

% Save the figure
saveas(gcf, '/home/ethan/fig_GLV(p3)/within_channel_syn_signifcoi_significant_channels_buttonpress_trials_GLV(p3).png');
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
sgtitle('Significant Redundant Co-I Values for Significant Channels Button presses as trials 100 Permutations', 'FontSize', 10);

% Save the figure
saveas(gcf, '/home/ethan/fig_GLV(p3)/within_channel_red_signifcoi_significant_channels_buttonpress_trials_GLV(p3).png');
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
sgtitle('Co-Information Heatmaps for Significant Channels Button presses as trials');

% Save the figure
saveas(gcf, '/home/ethan/fig_GLV(p3)/within_channel_coi_significant_channels_buttonpress_trials_GLV(p3).png');
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
%% number of trials identified 
disp('length of l_ vector created');
disp(length(l_));
disp(l_);
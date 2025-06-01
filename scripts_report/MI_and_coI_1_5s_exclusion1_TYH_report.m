clear 
clc
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
load('/home/ethan/Data/TYH_volitional_streaming_04.mat')

%% finding the indices where perception is 1 stream or 2 streams excluding all the values that are within 1.5s of each button press

% Defining the exclusion window in terms of index length
exclusionWindow = 1.5 * srate; % 1.5 seconds

% Find all indices for 1 and 2
indices_1 = find(perception == 1);
indices_2 = find(perception == 2);

% Force indices to row vectors
indices_1 = indices_1(:)';
indices_2 = indices_2(:)';

% Identify block boundaries for 1
diff_1 = [true, diff(indices_1) > 1, true];
block_starts_1 = find(diff_1(1:end-1));
block_ends_1 = find(diff_1(2:end)) - 1;

% Identify block boundaries for 2
diff_2 = [true, diff(indices_2) > 1, true];
block_starts_2 = find(diff_2(1:end-1));
block_ends_2 = find(diff_2(2:end)) - 1;

% Filter the indices for 1
filtered_indices_1 = [];
for i = 1:length(block_starts_1)
    block = indices_1(block_starts_1(i):block_ends_1(i));
    filtered_indices_1 = [filtered_indices_1, block((exclusionWindow + 1):end-exclusionWindow)];
end

% Filter the indices for 2
filtered_indices_2 = [];
for i = 1:length(block_starts_2)
    block = indices_2(block_starts_2(i):block_ends_2(i));
    filtered_indices_2 = [filtered_indices_2, block((exclusionWindow + 1):end-exclusionWindow)];
end

% Display results
disp(['Length of filtered_indices_1: ', num2str(length(filtered_indices_1))]);
disp(['Length of filtered_indices_2: ', num2str(length(filtered_indices_2))]);

%% find the corresponding values in stim for these indices in oneIndices and twoIndices

% Initialize vectors with 3s, same length as stim
stimOne = 3 * ones(size(stim)); 
stimTwo = 3 * ones(size(stim));

% found values in stim retain the same indices
stimOne(filtered_indices_1) = stim(filtered_indices_1); % Retain original indices
stimTwo(filtered_indices_2) = stim(filtered_indices_2); % Retain original indices

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
selectedIndices_stimOne = startIndices_stimOne(1:3:end); 
selectedIndices_stimTwo = startIndices_stimTwo(1:3:end); 

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
    [ii_, zii_] = ii_decode(l_, rr_, 't', n_t, 'pl', n_pl3, 'is_x', false);

    
    % Store results directly into the main cell arrays
    ii_all{ch} = ii_;
    zii_all{ch} = zii_;
end


%% 
% Prepare inputs for get_plotting_CoI_simple

% Preallocate arrays
%T = size(ii_all{find(~cellfun(@isempty, ii_all), 1)}, 1); % Timepoints
T = 721; % Number of timepoints
E = length(significant_channels); % Number of significant channels

% Initialize matrices
FFi = zeros(T, T, E);     % Raw CoI
FFm = zeros(T, T, E);     % Significance mask
FFmr = zeros(T, T, E);    % Redundancy mask
FFms = zeros(T, T, E);    % Synergy mask
MI1 = zeros(T, E);        % MI1
MI2 = zeros(T, E);        % MI2

% Fill the matrices
for i = 1:E
    ch = significant_channels(i);

    FFi(:, :, i) = ii_all{ch}.coi;
    FFm(:, :, i) = zii_all{ch}.coi;
    FFmr(:, :, i) = zii_all{ch}.red;
    FFms(:, :, i) = zii_all{ch}.syn;
    MI1(:, i) = ii_all{ch}.mi;
    MI2(:, i) = ii_all{ch}.mi; % assuming symmetric for plotting
end

% Pack into structure
acc.atDatasig.FFi = FFi;
acc.atDatasig.FFm = FFm;
acc.atDatasig.FFmr = FFmr;
acc.atDatasig.FFms = FFms;
acc.atDatasig.MI1 = MI1;
acc.atDatasig.MI2 = MI2;

% Now you can pass `acc` to get_plotting_CoI_simple
[FFIall, FFi_mean, FFi_r, FFi_s, FFm_sum, FFmr_sum, FFms_sum, MI1_out, MI2_out] = get_plotting_CoI_simple1(acc);


%% plotting 
% Axis and label parameters
xticks_CoI = [0, 275, 550];
xlimits = [-50 550];
x_labels = {'0', '275', '550'};
yticks_CoI = [0, 275, 550];
yticks_MI = [0, round(max(MI2_out(:))+0.25*max(MI2_out(:)), 3, 'significant') ];
ylimits_CoI = [-50 550];        % Adjust based on data range
ylimit_MI = [round(min(MI2_out(:)) + 0.1 * min(MI2(:)), 3, 'significant'), round(max(MI2(:)) + 0.25 * max(MI2(:)), 3, 'significant')];  % Adjust based on MI data range
y_labels_CoI = [0, 275, 550];
y_labels_MI = [0, round(max(MI2_out(:)) + 0.25 * max(MI2_out(:)), 3, 'significant')];
timing = (-60:660) * (1 / 1.2);  % gives you timing in milliseconds
temp_nb = length(significant_channels); % Number of significant channels/number im averaging over
%temp_nb = ['out of ' num2str(length(significant_channels))]; % Number of significant channels/number im averaging over
% Color limits
climits = [-fix(1000*max(abs(FFi_mean(:))))/1000, fix(1000*max(abs(FFi_mean(:))))/1000]; % Truncated to 3 decimal places
climits_mask = [0 1];         % Example; set based on your masking sets limits for signifcant plots

CoI_figure = plot_temporal(FFi_mean, FFm_sum, FFmr_sum, FFms_sum, MI1, MI2, FFi_r, FFi_s,  temp_nb, timing, ... 
            xlimits, ylimits_CoI, ylimit_MI, xticks_CoI, yticks_CoI, yticks_MI, x_labels, y_labels_CoI, y_labels_MI,climits, climits_mask);


% Add super title
sgtitle('TYH', 'FontSize', 25 );

% Save the figure
saveas(gcf, '/home/ethan/fig_report/MI_and_coI_1_5s_exclusion_plotted_together_TYH(p1).png');







%% Display results
disp('Significant channels:');
disp(significant_channels);
disp('number of channels with significant MI values');
disp(length(significant_channels));

%% percentage of triplets identified using exclusion method
Percentage_triplets_identified = length(combinedIndices) / 500 * 100;
disp('Percentage of triplets retained using direct extraction method');
disp(Percentage_triplets_identified);

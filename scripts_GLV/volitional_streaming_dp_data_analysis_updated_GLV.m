% load the data 

clear
load('/home/ethan/Data/GLV_volitional_streaming_04.mat')
whos
%% checking that each section of 0s in stim are the same length
% Find indices where stim is 0
zeroIndices = find(stim == 0);

% Identify breaks between separate sections
diffs = diff(zeroIndices);
breaks = find(diffs > 1);

% Ensure startIndices and endIndices are defined properly
if isempty(breaks)  % Case where all 0s are in one continuous block
    startIndices = zeroIndices(1);
    endIndices = zeroIndices(end);
else
    startIndices = [zeroIndices(1); zeroIndices(breaks + 1)];
    endIndices = [zeroIndices(breaks); zeroIndices(end)];
end

% Ensure that startIndices and endIndices are column vectors
startIndices = startIndices(:);
endIndices = endIndices(:);

% Compute lengths of each section
sectionLengths = endIndices - startIndices + 1;

% Check if all section lengths are the same
if all(sectionLengths == sectionLengths(1))
    disp('All sections of 0s have the same length.');
else
    disp('Sections of 0s have different lengths.');
    disp('Section lengths:');
    disp(sectionLengths);
    disp('Number of sections:');
    disp(length(sectionLengths));
end

%% checking that each section of 1s in stim are the same length
% Find indices where stim is 1
oneIndices = find(stim == 1);

% Identify breaks between separate sections
diffs = diff(oneIndices);
breaks = find(diffs > 1);

% Ensure startIndices and endIndices are defined properly
if isempty(breaks)  % Case where all 1s are in one continuous block
    startIndices = oneIndices(1);
    endIndices = oneIndices(end);
else
    startIndices = [oneIndices(1); oneIndices(breaks + 1)];  
    endIndices = [oneIndices(breaks); oneIndices(end)];
end

% Ensure that startIndices and endIndices are column vectors
startIndices = startIndices(:);
endIndices = endIndices(:);

% Compute lengths of each section
sectionLengths = endIndices - startIndices + 1;

% Check if all section lengths are the same
if all(sectionLengths == sectionLengths(1))
    disp('All sections of 1s have the same length.');
else
    disp('Sections of 1s have different lengths.');
    disp('Section lengths:');
    disp(sectionLengths);
    disp('number of sections:');
    disp(length(sectionLengths));
end
uniqueValues = unique(stim);
disp('Unique values in stim:');
disp(uniqueValues);
%% checking that each section of 1s in stim are the same length
% Find indices where stim is 1
oneIndices = find(perception == 1);

% Identify breaks between separate sections
diffs = diff(oneIndices);
breaks = find(diffs > 1);

% Ensure startIndices and endIndices are defined properly
if isempty(breaks)  % Case where all 1s are in one continuous block
    startIndices = oneIndices(1);
    endIndices = oneIndices(end);
else
    startIndices = [oneIndices(1); oneIndices(breaks + 1)];  
    endIndices = [oneIndices(breaks); oneIndices(end)];
end

% Ensure that startIndices and endIndices are column vectors
startIndices = startIndices(:);
endIndices = endIndices(:);

% Compute lengths of each section
sectionLengths = endIndices - startIndices + 1;

% Check if all section lengths are the same
if all(sectionLengths == sectionLengths(1))
    disp('All sections of 1s have the same length.');
else
    disp('Sections of 1s have different lengths.');
    disp('Section lengths:');
    disp(sectionLengths);
    disp('number of 1 sections:');
    disp(length(sectionLengths));
end
%% checking that each section of 1s in stim are the same length
% Find indices where stim is 1
oneIndices = find(perception == 2);

% Identify breaks between separate sections
diffs = diff(oneIndices);
breaks = find(diffs > 1);

% Ensure startIndices and endIndices are defined properly
if isempty(breaks)  % Case where all 1s are in one continuous block
    startIndices = oneIndices(1);
    endIndices = oneIndices(end);
else
    startIndices = [oneIndices(1); oneIndices(breaks + 1)];  
    endIndices = [oneIndices(breaks); oneIndices(end)];
end

% Ensure that startIndices and endIndices are column vectors
startIndices = startIndices(:);
endIndices = endIndices(:);

% Compute lengths of each section
sectionLengths = endIndices - startIndices + 1;

% Check if all section lengths are the same
if all(sectionLengths == sectionLengths(1))
    disp('All sections of 2s have the same length.');
else
    disp('Sections of 2s have different lengths.');
    disp('Section lengths:');
    disp(sectionLengths);
    disp('number of 2 sections:');
    disp(length(sectionLengths));
end
uniqueValues = unique(perception);
disp('Unique values in perception:');
disp(uniqueValues);
%% sections are found to be of lengths 145, 144, need to make them
% all length 140

targetLength = 145;
endIndices(sectionLengths == 144) = startIndices(sectionLengths == 144) + targetLength - 1;
endIndices(sectionLengths == 142) = startIndices(sectionLengths == 142) + targetLength - 1;
% Compute new section lengths after modification
newSectionLengths = endIndices - startIndices + 1;

% Verify all sections are exactly 140
assert(all(newSectionLengths == targetLength), 'Error: Not all sections are of length 140.');

% Print confirmation
disp('✅ All sections have been successfully truncated to 145 rows.');
 %disp('Final section lengths:');
 %disp(newSectionLengths);

%% Filter data using EEGlab

% Channels on first dimension
dp_data = dp_data'; 

% import data to EEGlab format
EEG = pop_importdata('dataformat','matlab','nbchan',0,'data',dp_data,'setname','data','srate',1200,'pnts',0,'xmin',0);

% Filter data 0.5-40 Hz
EEG = pop_eegfiltnew(EEG, 'locutoff',0.5,'hicutoff',40,'plotfreqz',1);

% back to channels second dimension
dp_data = EEG.data';

%% transforming the 2D dp_data matrix into a 3d matrix with the 3rd dim being each section of 1s found in stim

% Define the number of columns in the original 2D matrix
numCols = size(dp_data, 2);

% Number of identified sections
numSections = length(startIndices);

% Preallocate a 3D matrix (140 rows, numCols columns, numSections depth)
dp_data_3D = zeros(140, numCols, numSections);

% Extract sections and store in 3D matrix
for i = 1:numSections
    dp_data_3D(:,:,i) = dp_data(startIndices(i):endIndices(i), :);
end

% Display confirmation
disp('✅ Successfully transformed the 2D matrix into a 3D matrix.');
disp(['Size of 3D matrix: ', num2str(size(dp_data_3D,1)), ' x ', num2str(size(dp_data_3D,2)), ' x ', num2str(size(dp_data_3D,3))]);

%% averaging over the 3rd dim to obtain a 2D matrix again

% Compute the mean over the 3rd dimension (averaging across sections)
dp_data_3D_avrg = mean(dp_data_3D, 3);

% Ensure the result is truly 2D by removing the singleton dimension
dp_data_3D_avrg = squeeze(dp_data_3D_avrg);

% Display confirmation
disp('✅ Successfully collapsed the 3D matrix into a true 2D matrix by averaging.');
disp(['Size of new 2D matrix: ', num2str(size(dp_data_3D_avrg,1)), ' x ', num2str(size(dp_data_3D_avrg,2))]);

%% plotting acouple of the channels (columns) as a sanity check

% Choose which column to plot 
columnIndex = 175;

% Extract the selected column
columnData = dp_data_3D_avrg(:, columnIndex);

% Create the x-axis (assuming time or row index)
x = 1:length(columnData);
time = x/srate;
% Plot the selected column
figure; % Create a new figure window
plot(time, columnData, 'b-', 'LineWidth', 2); % Blue line plot

% Add labels and title
xlabel('time (s)'); % X-axis label
ylabel('voltage '); % Y-axis label
title(['Plot of Channel ', num2str(columnIndex)]); 

%% plot of all channels for comparison without error bars
% Get the number of channels
numChannels = size(dp_data_3D_avrg, 2);

% Define grid layout
gridRows = ceil(sqrt(numChannels)); 
gridCols = ceil(numChannels / gridRows);

% Define x-axis (assuming time or samples)
x = 1:size(dp_data_3D_avrg, 1);
time = x/srate;
% Create figure
figure;

% Loop through each channel and create subplots in a grid layout
for i = 1:numChannels
    subplot(gridRows, gridCols, i); % Arrange in a grid
    plot(time, dp_data_3D_avrg(:, i), 'b', 'LineWidth', 1); % Plot channel in blue
    
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
    yLimit = max(abs(dp_data_3D_avrg(:))); % Get max absolute value for symmetry
    ylim([-yLimit, yLimit]); % Make y-axis symmetric around 0
    
    % Add smaller channel number label in the top-right corner
    text(max(time) * 0.95, max(dp_data_3D_avrg(:)) * 0.9, num2str(i), ...
        'FontSize', 7, 'FontWeight', 'Bold', ...
        'Units', 'data', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top');
    
    
end

% Add super title
sgtitle('ERP of each Channel'); 

%% plot with SEM (error bars, standard error of the mean)

% Channels in first dimension
dp_data_ep = permute(dp_data_3D, [2 1 3]);

figure;

for chan = 1:size(dp_data_ep,1)
subplot(10,20,chan)
 stdshade_acj(squeeze(dp_data_ep(chan,:,:))',0.2,'r')

end

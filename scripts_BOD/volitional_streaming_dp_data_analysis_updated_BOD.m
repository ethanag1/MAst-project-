clear
clc
% load the data 
load('/home/ethan/Data/BOD_volitional_streaming_04.mat')
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
%% sections are found to be of lengths 140, 141 and 142, need to make them
% all length 140

targetLength = 140;
endIndices(sectionLengths == 141) = startIndices(sectionLengths == 141) + targetLength - 1;
endIndices(sectionLengths == 142) = startIndices(sectionLengths == 142) + targetLength - 1;
% Compute new section lengths after modification
newSectionLengths = endIndices - startIndices + 1;

% Verify all sections are exactly 140
assert(all(newSectionLengths == targetLength), 'Error: Not all sections are of length 140.');

% Print confirmation
disp('✅ All sections have been successfully truncated to 140 rows.');
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



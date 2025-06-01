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

%% histogram plot of each 

% Define the number of bins
numBins = 20;

% Determine common bin edges based on the combined range
combinedMin = min([IRI_1stream; IRI_2stream]);
combinedMax = max([IRI_1stream; IRI_2stream]);
binEdges = linspace(combinedMin, combinedMax, numBins+1); % Ensure same bin size for both


% Create a figure
figure;

% First histogram (Red)
subplot(1,2,1); % First subplot
histogram(IRI_1stream, 'BinEdges', binEdges, 'FaceColor', 'r'); % Use predefined bin edges
title('IRI One-Stream');
xlabel('Inter-response interval (S)');
ylabel('Occurrences');
hold on;

% Second histogram (Blue)
subplot(1,2,2); % Second subplot
histogram(IRI_2stream, 'BinEdges', binEdges, 'FaceColor', 'b'); % Use same bin edges
title('IRI Two-Stream');
xlabel('Inter-response interval (S)');
ylabel('Occurrences');
hold on;

% Set same x and y limits
subplot(1,2,1);
xlim([combinedMin combinedMax]); % Set consistent x-axis limits
ylim([0 max([histcounts(IRI_1stream, binEdges), histcounts(IRI_2stream, binEdges)], [], 'all')]); % Consistent y-axis

subplot(1,2,2);
xlim([combinedMin combinedMax]); % Set consistent x-axis limits
ylim([0 max([histcounts(IRI_1stream, binEdges), histcounts(IRI_2stream, binEdges)], [], 'all')]); % Consistent y-axis

% Add text with statistics in top right, displaying Count first
subplot(1,2,1);
text(max(xlim)-1, max(ylim)-1, sprintf('N: %d\nMean: %.2f\nMedian: %.2f', N1, meanOne, medianOne), ...
    'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'FontSize', 10, 'FontWeight', 'bold');

subplot(1,2,2);
text(max(xlim)-1, max(ylim)-1, sprintf('N: %d\nMean: %.2f\nMedian: %.2f', N2, meanTwo, medianTwo), ...
    'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'FontSize', 10, 'FontWeight', 'bold');

% Save the figure
saveas(gcf, '/home/ethan/fig_GLV(p3)/IRI_histograms_1_and_2_seperate_GLV(p3).png');

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

%% Histogram plot of combined IRI

% Plot histogram
figure;
histogram(combinedLengths, 20, 'FaceColor', 'g', 'EdgeColor', 'k'); % Green bars, black edges
title('Combined IRI');
xlabel('Inter-response interval (S)');
ylabel('Occurrences');

% Add text with statistics in the top-right
xLimits = xlim;
yLimits = ylim;
text(xLimits(2) - 1, yLimits(2) - 1, ...
    sprintf('Total: %d\nMean: %.2f\nMedian: %.2f', numCombined, meanCombined, medianCombined), ...
    'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', ...
    'FontSize', 10, 'FontWeight', 'bold', 'Color', 'k');

    % Save the figure
saveas(gcf, '/home/ethan/fig_GLV(p3)/IRI_histogram_1_and_2_Combined_GLV(p3).png');
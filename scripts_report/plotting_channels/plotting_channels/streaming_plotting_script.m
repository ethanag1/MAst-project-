
clear 
clc
opengl('save', 'software')
opengl software
pathToBrain = '/home/ethan/Data/BOD_brain.mat'; %enter path to brain
addpath(genpath('/home/ethan/scripts_report/spm_25.01.02/spm'))

load(pathToBrain)

figure
 
surf = ieeg_RenderGifti(cortex); % eneter cortex_R/L or cortex to plot hemispheres of choice
surf.FaceAlpha = 0.4; % set transparency of the brain %0.2

ieeg_viewLight(90, 0); % change the angle of view  -100, 40 -90,0 leftview , 90,0 right view , 0,90 top view

% since we want to plot all electrodes as circles of the same diameter we assign a weight of 1 to all
wts = ones(size(dp_locs,1),1); % change "dp_locs" to "locs" if you want raw electrode positions and not bipolar channels

plotCortexWeights(dp_locs, wts) % change "dp_locs" to "locs" if you want raw electrode positions and not bipolar channels
% open this function to change the diameter (line 6, msize) and color of circles (lines 24-30)

%saveas(figure, '/home/ethan/fig_report/electrode_brain_mapping_TYH.png');
exportgraphics(gcf, '/home/ethan/fig_report/electrode_brain_mapping_offcentre_BOD1.png', 'Resolution', 600);

disp('figure created');
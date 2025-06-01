function [FFIall, FFi_mean, FFi_r, FFi_s, FFm_sum, FFmr_sum, FFms_sum, MI1, MI2] = get_plotting_CoI_simple1(acc)
%GET_PLOTTING_COI_SIMPLE Simplified CoI plotting data extractor.
%   This function replaces get_plotting_CoI_ALL for cases where all
%   electrodes come from the same area. It takes a single structure
%   `acc` with field `atDatasig` and computes summary CoI measures.
%
%   Inputs:
%       acc: struct containing .atDatasig with fields:
%           - MI1    : [T x E] mutual information channel 1
%           - MI2    : [T x E] mutual information channel 2
%           - FFi    : [T x T x E] raw CoI matrices per electrode/pair
%           - FFm    : [T x T x E] binary significance masks
%           - FFmr   : [T x T x E] redundant-significance masks
%           - FFms   : [T x T x E] synergetic-significance masks
%
%   Outputs:
%       FFIall   : [T x T x E] raw CoI stacks (same as acc.atDatasig.FFi)
%       FFi_mean : [T x T] mean CoI across electrodes
%       FFi_r    : [T x T] redundant component (>=0) of mean CoI
%       FFi_s    : [T x T] synergetic component (<=0) of mean CoI
%       FFm_sum  : [T x T] sum of all significance masks
%       FFmr_sum : [T x T] sum of redundant masks
%       FFms_sum : [T x T] sum of synergetic masks
%       MI1      : [T x E] mutual info channel 1 (from acc.atDatasig.MI1)
%       MI2      : [T x E] mutual info channel 2 (from acc.atDatasig.MI2)

% Extract raw stacks and MI vectors
FFIall = acc.atDatasig.FFi;
MI1    = acc.atDatasig.MI1;
MI2    = acc.atDatasig.MI2;

% Compute mean CoI across the third dimension (electrodes/pairs)
FFi_mean = squeeze(mean(FFIall, 3));

% Split into redundant (non-negative) vs. synergetic (non-positive)
FFi_r = FFi_mean;
FFi_r(FFi_mean <= 0) = 0;
FFi_s = FFi_mean;
FFi_s(FFi_mean >= 0) = 0;

% Sum up significance masks across electrodes/pairs
FFm_sum  = sum(logical(acc.atDatasig.FFm),3);
FFmr_sum = sum(logical(acc.atDatasig.FFmr),3);
FFms_sum = sum(logical(acc.atDatasig.FFms),3);
end

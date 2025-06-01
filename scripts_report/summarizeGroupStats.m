function [FFi_mean_group, FFi_r_group, FFi_s_group, ...
          FFm_sum_group, FFmr_sum_group, FFms_sum_group, ...
          MI1_group, MI2_group] = summarizeGroupStats( ...
          groupFFi, groupFFm, groupFFmr, groupFFms, groupMI1, groupMI2)
% SUMMARIZEGROUPSTATS  Compute group-level summaries of FF and MI data
    % mean of FF indexing across 3rd dimension
    FFi_mean_group = mean(groupFFi, 3);
    % positive and negative parts
    FFi_r_group    = max(FFi_mean_group, 0);
    FFi_s_group    = min(FFi_mean_group, 0);

    % sum of FF magnitude masks across subjects
    FFm_sum_group  = sum(logical(groupFFm),  3);
    FFmr_sum_group = sum(logical(groupFFmr), 3);
    FFms_sum_group = sum(logical(groupFFms), 3);

    % MI time series remain unchanged
    MI1_group = groupMI1;
    MI2_group = groupMI2;
end

%Plots temporal electrodes 
%saves figures as Pdf and fig
function [CoI_figure] = plot_temporal(tempFFi, tempFFm, tempFFmr, tempFFms, tempFFmi1, tempFFmi2, tempFFir, tempFFis, temp_nb, timing, ...
                                        xlimits,ylimits_CoI,ylimit_MI,xticks_CoI,yticks_CoI,yticks_MI,x_labels,y_labels_CoI, y_labels_MI,climits,climits_mask)

%% Create figure and plot
CoI_figure = figure(1);
tiledlayout(4,3,'TileSpacing','Compact');

% This is the combined synergetic/redundant
ax(1) = nexttile(5);
set(gcf,'renderer','Painters','Position', [10 10 800 450]);
contourf(timing,timing,tempFFi, 50, 'LineStyle', 'none');
hold on
xlim(xlimits); ylim(ylimits_CoI); clim(climits)
set(gca,'ytick',yticks_CoI, 'yticklabel', yticks_CoI, 'xtick',xticks_CoI, 'xticklabel', xticks_CoI);
colormap(ax(1), redblue)
set(gca, 'box', 'on'); % Restore border
set(gca, 'XColor', 'k', 'YColor', 'k', 'LineWidth', 1.2) % Restore tick marks and border color

% MI1
ax(12) = nexttile(2);
stdshade_acj(tempFFmi1',.2,'k', timing);
xlim(xlimits); ylim(ylimit_MI);
set(gca,'ytick',yticks_MI, 'yticklabel', yticks_MI, 'xtick', xticks_CoI, 'xticklabel', xticks_CoI);
set(gca, 'box', 'off');
ylabel('MI (bits)', 'Rotation', 90)
xlabel('Time (ms)') % NEW: X-axis label

% MI2
ax(14) = nexttile(4);
stdshade_acj(tempFFmi2',.2,'k',timing); view(-90,90)
xlim(xlimits); ylim(ylimit_MI);
set(gca,'ytick',yticks_MI, 'yticklabel', y_labels_MI, 'xtick', xticks_CoI, 'xticklabel', x_labels);
set(gca, 'XAxisLocation', 'top');
set(gca, 'box', 'off');
ylabel('MI (bits)')
xlabel('Time (ms)') % NEW: Vertical Y-axis label

% Combined synergetic/redundant mask
ax(2) = nexttile(6);
imagesc(timing,timing,tempFFm./temp_nb);
title(num2str(temp_nb));
set(gca,'YDir','normal');
xlim(xlimits); ylim(ylimits_CoI); clim(climits)
set(gca,'ytick',yticks_CoI, 'yticklabel', y_labels_CoI, ...
        'xtick',xticks_CoI, 'xticklabel', x_labels, ...
        'clim', climits_mask);
colormap(ax(2), flipud(bone))
xlabel('Time (ms)') % NEW: X-axis label
ylabel('Time (ms)') % NEW: Y-axis label

% Redundancy CoI
ax(3) = nexttile(8);
set(gcf,'renderer','Painters','Position', [10 10 800 450]);
contourf(timing,timing,tempFFir, 50,'linecolor','none');
xlim(xlimits); ylim(ylimits_CoI); clim(climits)
set(gca,'ytick',yticks_CoI, 'yticklabel', y_labels_CoI, ...
        'xtick',xticks_CoI, 'xticklabel', x_labels);
colormap(ax(3), redblue)
ylabel({'Redundancy', 'only'}, 'Color', 'r', 'Rotation', 90, 'FontSize', 10);

% Synergy CoI
ax(4) = nexttile(11);
set(gcf,'renderer','Painters','Position', [10 10 800 450]);
contourf(timing,timing,tempFFis, 50,'linecolor','none');
xlim(xlimits); ylim(ylimits_CoI); clim(climits)
set(gca,'ytick',yticks_CoI, 'yticklabel', y_labels_CoI, ...
        'xtick',xticks_CoI, 'xticklabel', x_labels);
colormap(ax(4), redblue)
ylabel('Synergy only', 'Color', 'b', 'Rotation', 90, 'FontSize', 10);
xlabel('Time (ms)') % NEW: X-axis label

% Redundant mask
ax(6) = nexttile(9);
imagesc(timing,timing,tempFFmr./temp_nb);
set(gca,'YDir','normal');
xlim(xlimits); ylim(ylimits_CoI);
set(gca,'ytick',yticks_CoI, 'yticklabel', y_labels_CoI, ...
        'xtick',xticks_CoI, 'xticklabel', x_labels, ...
        'clim', climits_mask);
colormap(ax(6), flipud(bone))
cb1 = colorbar;
cb1.Ticks = [0 1]; % Explicit ticks
ylabel(cb1, 'Fraction', 'FontSize', 10);
ylabel('Time (ms)') % NEW: Y-axis label

% Synergy mask
ax(7) = nexttile(12);
imagesc(timing,timing,tempFFms./temp_nb);
set(gca,'YDir','normal');
xlim(xlimits); ylim(ylimits_CoI);
set(gca,'ytick',yticks_CoI, 'yticklabel', y_labels_CoI, ...
        'xtick',xticks_CoI, 'xticklabel', x_labels, ...
        'clim', climits_mask);
colormap(ax(7), flipud(bone))
xlabel('Time (ms)') % NEW: X-axis label
ylabel('Time (ms)') % NEW: Y-axis label

% Colorbar Co-I (Tile 1)
ax(8) = nexttile(1);
set(gca,'ytick',[], 'yticklabel', [], 'xtick',[], 'xticklabel', [], ...
    'clim', climits, 'box', 'off', ...
    'XColor', 'none', 'YColor', 'none', 'LineWidth', 0.01); % Fully remove visible borders
cb2 = colorbar;
cb2.Ticks = [climits(1) 0 climits(2)]; % Explicit ticks
ylabel(cb2, 'co-I (bits)', 'Rotation', 90, 'FontSize', 10);
colormap(ax(8), redblue);
hold off

% Tick mark formatting and border thickness
for i = 1:length(ax)
    if isgraphics(ax(i))
        set(ax(i), 'TickDir', 'out', 'TickLength', [0.02 0.02], ...
            'LineWidth', 1.2)
    end
end

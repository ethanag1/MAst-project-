function plotCortexWeights(locs, wts, threshFrac)
%
% Modified by kjm & HH 2021
%

    msize=3;
    if nargin < 4, threshFrac = 0.1; end % abs(values) below threshFrac*wm deemed not significant

    % wts(isnan(wts)) = 0;
    
    % this function plots colored dots in brain rendering
    wm=max(abs(wts));
    pthresh = threshFrac*wm;

    hold on
    for k=1:size(locs,1)
        if abs(wts(k)) == 0 % insignificant
            plot3((locs(k,1)),locs(k,2),locs(k,3),'o',... % plot3(-abs(locs(k,1)),locs(k,2),locs(k,3),'o',...to plot left Hemi only 
                'MarkerSize',msize*abs(wts(k))/wm+msize,...
                'LineWidth',0.1,...
                'MarkerEdgeColor',[0.5 0.5 0.5],... 
                'MarkerFaceColor',[0 0 1])  
            
        elseif wts(k) == 1 % positively significant
 
            plot3((locs(k,1)),locs(k,2),locs(k,3),'o',...
                'MarkerSize',msize*abs(wts(k))/wm+msize,...
                'LineWidth',0.1,...
                'MarkerEdgeColor',[0.5 0.5 0.5],...
                'MarkerFaceColor',[1 0 0])  
            
            %text((locs(k,1)),locs(k,2),locs(k,3),pts{k})

        elseif wts(k) == 2 % positively significant

            plot3((locs(k,1)),locs(k,2),locs(k,3),'ro',...
                'MarkerSize',msize,... 
                'LineWidth',0.1,...
                'MarkerEdgeColor',[0.5 0.5 0.5],...
                'MarkerFaceColor',[0 1 1]) 
        end
    end
    hold off
    
end

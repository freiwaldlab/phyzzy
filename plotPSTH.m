function [ psthAxes ] = plotPSTH(psthArray, psthAxes, psthParams, plotType, psthTitle, ylabels)
%UNTITLED10 Summary of this function goes here
%   Detailed explanation goes here

 
psthPre = psthParams.psthPre;
psthImDur = psthParams.psthImDur;
psthPost = psthParams.psthPost;

xrange= [-psthPre psthImDur+psthPost]; 
nrows = size(psthArray,1);
yaxis = [1 nrows];
title(psthTitle); 
if strcmp(plotType,'color')
  try
    caxis([min(min(psthArray)),max(max(psthArray))]);
  catch
    disp('start error message');
    disp(psthArray);
    disp([min(min(psthArray)),max(max(psthArray))]);
    assert(0,'failed in plotPSTH');
    return
  end
  if isempty(psthAxes)
    psthAxes = imagesc(xrange, yaxis, psthArray);
  else
    imagesc(psthAxes, xrange, yaxis, psthArray);
  end
  h = colorbar;
  ylabel(h,'Firing Rate [Hz]','FontSize',14);
  colormap(psthParams.colormap);
  ylimits= ylim();
  yRange = ylimits(2) - ylimits(1);
  hold on
  if (psthPre+psthPost)/psthImDur > 20
    stimDurLineWidth = 0.1;
  else
    stimDurLineWidth = 4;
  end
  vertLineColor = [0.5, 0.5, 0.5];
  draw_vert_line(0,'Color',vertLineColor,'LineWidth',stimDurLineWidth);
  draw_vert_line(psthImDur,'Color',vertLineColor,'LineWidth',stimDurLineWidth);
  set(gca,'YTick',linspace(ylimits(1)+yRange/(2*nrows),ylimits(2)-yRange/(2*nrows),nrows),'YTicklabel',ylabels,...
    'box','off','TickDir','out','FontSize',14,'TickLength',[.012 .012],'LineWidth',1.3);
else
  error('psth line plot not yet implemented');
end
xlabel('Time from stimulus onset (ms)', 'FontSize',14);
hold off

Output.DEBUG('min of psth %s: %d\n',psthTitle,(min(min(psthArray))));
Output.DEBUG('max of psth %s: %d\n',psthTitle,(max(max(psthArray))));
end


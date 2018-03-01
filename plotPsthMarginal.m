function [ psthAxes ] = plotPsthMarginal( psthArray, psthPre, psthPost, psthImDur, frCalcWin, psthTitle, ylabels, psthColormap )
%UNTITLED10 Summary of this function goes here
%   Detailed explanation goes here

xrange= [-psthPre psthImDur+psthPost]; 
nrows = size(psthArray,1);
yrange= [1 nrows];

ROWS = 4;
PSTH_COLS = 12; % cols is short for columns
MARGINAL_COLS = 4;
%cols is the number of subplot columns; the +2 is for the extra colorbar
%column, and spacer between the color plot and the marginals
cols = PSTH_COLS+MARGINAL_COLS+2; 

try
  caxis([min(min(psthArray)),max(max(psthArray))]);
catch
  disp('start error message');
  disp(psthArray);
  disp([min(min(psthArray)),max(max(psthArray))]);
  assert(0,'failed in plotPSTH');
  return
end

% calculate the marginals (todo: give option for passed in, with errorbars)
marginals = sum(psthArray(:,psthPre+1+frCalcWin(1):psthPre+1+frCalcWin(2)),2);

% set up the subplots
psthSlots = zeros(ROWS,PSTH_COLS);
for row_i = 1:ROWS
  psthSlots(row_i,:) = ((row_i - 1)*cols+1):((row_i - 1)*cols+PSTH_COLS);
end
psthSlots = reshape(psthSlots,numel(psthSlots),1);

marginalSlots = zeros(ROWS,MARGINAL_COLS);
for row_i = 1:ROWS
  marginalSlots(row_i,:) = ((row_i - 1)*cols+PSTH_COLS+2):((row_i - 1)*cols+PSTH_COLS+1+MARGINAL_COLS);
end
marginalSlots = reshape(marginalSlots,numel(marginalSlots),1);

colorbarSlots = cols*(1:ROWS);

figure();
% Plot the heatmap psth
ah1 = subplot(ROWS,cols,psthSlots); 

psthAxes = imagesc(xrange, yrange, psthArray);
colormap(psthColormap);
ylimits= ylim();
yRange = ylimits(2) - ylimits(1);  %todo: remove egregious use of this variable twice, but with different caps
hold on
if (psthPre+psthPost)/psthImDur > 20
  stimDurLineWidth = 0.1;
else
  stimDurLineWidth = 4;
end
draw_vert_line(0,'Color',[0.8,0.8,0.9],'LineWidth',stimDurLineWidth);
draw_vert_line(psthImDur,'Color',[0.8,0.8,0.9],'LineWidth',stimDurLineWidth);
set(gca,'YTick',linspace(ylimits(1)+yRange/(2*nrows),ylimits(2)-yRange/(2*nrows),nrows),'YTicklabel',ylabels,...
  'box','off','TickDir','out','FontSize',14,'TickLength',[.012 .012],'LineWidth',1.3);

xlabel('Time from stimulus onset (ms)', 'FontSize',14);
hold off

% plot the colorbar on the right edge of the graph
map = colormap(gca);
subplot(ROWS,cols,colorbarSlots);
colormap(map)
%plot([1,2,3,4,5],'linestyle','none');
set(gca,'Visible','off')
h = colorbar('location','WestOutside','AxisLocation','in');
set(gca,'FontSize',14);
ylabel(h,'Hz','FontSize',14);

% plot the marginals
ah2 = subplot(ROWS,cols,marginalSlots); %todo: set size to exactly match psth
plot(marginals,'linewidth',4,'marker','s','color','k');
set(gca,'CameraUpVector',[1,0,0])
set(gca,'YDir','reverse');
set(gca,'XTickLabel','');
set(gca,'FontSize',14);
ylabel('fr','FontSize',14)
% align y axes in rotated plot; see http://stackoverflow.com/questions/5259853/matlab-how-to-align-the-axes-of-subplots-when-one-of-them-contains-a-colorbar
% note: drawnow line synchronizes the rendering thread with the main
drawnow;
pos1 = get(ah1,'Position');
pos2 = get(ah2,'Position');
pos2(2) = max(pos1(2),pos2(2)); % right limit
pos1(2) = max(pos1(2),pos2(2));
pos2(4) = min(pos1(4),pos2(4)); % left limit
pos1(4) = min(pos1(4),pos2(4));
set(ah2,'Position',pos2);
set(ah1,'Position',pos1);
%%%
suptitle(psthTitle);

Output.DEBUG('min of psth %s: %d\n',psthTitle,(min(min(psthArray))));
Output.DEBUG('max of psth %s: %d\n',psthTitle,(max(max(psthArray))));
end


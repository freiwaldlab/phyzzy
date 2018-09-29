function [  ] = subplotsFromFigs( dir, fname, figures, equalColWidths, equalRowHeights )
%texSubplot uses the beamer class of LaTeX to make a subplot figure
%   If same fig on sequential rows, fig spans those rows
%   Params: 
%   - dir, string, location at which to store tmp pngs
%   - fname: full path+name of output png
%   - figures: cell array of cell arrays of figure handles, .fig filenames, or .png
%              filenames. The first cell array index gives the row index in the assembled figure; 
%               the second gives the column in the assembled final figure
%               When the same filename or handle appears in adjacent rows
%               and columns, the corresponding figure will span those rows
%               and columns
%   - equalColWidths, equalRowHeights: if 1, make all cols/rows the same size. If 0, cols/rows will
%               be sized according the width of their largest constituent figure


interRowSpace = 80;
interColSpace = 80;
verticalMargin = 120;
horizontalMargin = 120;
verticalSizes = zeros(length(figures),length(figures{1}));
horizontalSizes = zeros(length(figures),length(figures{1}));
rowsSpanned = ones(length(figures),length(figures{1}));
colsSpanned = ones(length(figures),length(figures{1}));
finished = zeros(length(figures),length(figures{1}));
figurePixels = cell(length(figures),length(figures{1}));
% preprocess the figure inputs, if not already png files
for fig_i = 1:length(figures)
  for fig_j = 1:length(figures{1})
    if finished(fig_i,fig_j)
      continue
    end
    % if a subplot is .fig or fig handle, generate png and update figures with its filename
    if ~ischar(figures{fig_i}{fig_j})  %it's a figure handle
      assert(isgraphics(figures{fig_i}{fig_j},'Figure'),'Invalid input type: choices are filename or figure handle');
      subplotFilename = sprintf('%s/subplot_%d_%d.png',dir,fig_i,fig_j);
      figure(figures{fig_i}{fig_j});
      set(gcf, 'color', [1 1 1]);
      export_fig(subplotFilename,'-m1.2','-opengl');
    elseif ~isempty(regexp(figures{fig_i}{fig_j},'.fig','ONCE'))  %it's a .fig
      h = openfig(figures{fig_i}{fig_j},'invisible');
      set(gcf, 'color', [1 1 1]);
      subplotFilename = regexprep(figures{fig_i}{fig_j},'.fig','.png');
      export_fig(subplotFilename,'-m1.2','-opengl');
      close(h);
    else
      subplotFilename = figures{fig_i}{fig_j};
    end
    pixels = imread(subplotFilename);
    figurePixels{fig_i}{fig_j} = pixels;
    verticalSizes(fig_i,fig_j) = size(pixels,1);
    horizontalSizes(fig_i,fig_j) = size(pixels,2);
    finished(fig_i,fig_j) = 1;
    % now, replace all the entries for this figure with the filename
    movesRight = 1;
    while fig_j+movesRight <= length(figures{1})
      if all(figures{fig_i}{fig_j+movesRight} == figures{fig_i}{fig_j}) %note: 'all' makes it work for strings and handles
        figures{fig_i}{fig_j+movesRight} = subplotFilename;
        finished(fig_i,fig_j+movesRight) = 1;
        movesRight = movesRight + 1;
      else
        break
      end
    end
    movesDown = 1;
    while fig_i+movesDown <= length(figures)
      if all(figures{fig_i+movesDown}{fig_j} == figures{fig_i}{fig_j})
        for moveRight = 0:movesRight-1
          finished(fig_i+movesDown,fig_j+moveRight) = 1;
          figures{fig_i+movesDown}{fig_j+moveRight} = subplotFilename;
        end
        movesDown = movesDown + 1;
      else
        break
      end
    end
    figures{fig_i}{fig_j} = subplotFilename;
    rowsSpanned(fig_i,fig_j) = movesDown;
    colsSpanned(fig_i,fig_j) = movesRight;
  end
end

if equalRowHeights
  rowHeights = ceil(max(max(verticalSizes./rowsSpanned - (rowsSpanned-1)*interRowSpace))*ones(size(verticalSizes,1),1));
else
  rowHeights = ceil(max(verticalSizes./rowsSpanned - (rowsSpanned-1)*interRowSpace,[],2));
end
if equalColWidths
  colWidths = ceil(max(max(horizontalSizes./colsSpanned - (colsSpanned-1)*interColSpace))*ones(size(verticalSizes,2),1));
else
  colWidths = ceil(max(horizontalSizes./colsSpanned - (colsSpanned-1)*interColSpace,[],1));
end

output = uint8(255*ones(sum(rowHeights)+interRowSpace*(length(rowHeights)-1)+2*verticalMargin,...
  sum(colWidths)+interColSpace*length(colWidths-1)+2*horizontalMargin,3));
for row_i = 1:length(figures)
  for col_i = 1:length(figures{1})
    if ~isempty(figurePixels{row_i}{col_i})
      pix = figurePixels{row_i}{col_i};
      topEdge = sum(rowHeights(1:row_i-1))+verticalMargin+(row_i-1)*interRowSpace;
      rightEdge = sum(colWidths(1:col_i-1))+horizontalMargin+(col_i-1)*interColSpace;
      output(topEdge:topEdge+size(pix,1)-1,rightEdge:rightEdge+size(pix,2)-1,:) = pix;
    end
  end
end
imwrite(output, fname);
end








function subplotsFromFigs(dir, fname, figures, params)
%subplotsFromFigs uses the beamer class of LaTeX to make a subplot figure
%   If same fig on sequential rows, fig spans those rows
%   Params:
%   - dir, string, location at which to store tmp pngs
%   - fname: full path+name of output png
%   - figures: cell array of cell arrays of figure handles, .fig filenames, or .png
%              filenames. The first cell array index gives the row index in the assembled figure;
%               the second gives the column in the assembled final figure
%               When the same filename or handle appears in adjacent rows
%               and columns, the corresponding figure will span those rows
%               and columns, e.g. {{fh1,fh2,fh3},{fh1,fh4,fh5}}, where the fhNs are any combination of figure handles and filenames
%   - equalColWidths, equalRowHeights: if 1, make all cols/rows the same size. If 0, cols/rows will
%               be sized according the width of their largest constituent figure

%% Default Variables
[equalColWidths, equalRowHeights] = deal(1);
interRowSpace = 10;
interColSpace = 10;
verticalMargin = 60;
horizontalMargin = 60;
pageColor = [1 1 1];

% Load variables from param if available
fields = fieldnames(params);
for param_i = 1:length(fields)
  eval(sprintf('%s = params.%s;', fields{param_i}, fields{param_i}));
end

%% Initialize needed variables
[verticalSizes, horizontalSizes, finished] = deal(zeros(length(figures),length(figures{1})));
[rowsSpanned, colsSpanned] = deal(ones(length(figures),length(figures{1})));
figurePixels = cell(length(figures),1);
[figurePixels{:}] = deal(cell(length(figures{1}),1));

for fig_i = 1:length(figures)
  for fig_j = 1:length(figures{1})
    if finished(fig_i,fig_j) || isempty(figures{fig_i}{fig_j})
      continue
    end
    % if a subplot is .fig or fig handle, generate png and update figures with its filename
    if ~ischar(figures{fig_i}{fig_j})  %it's a figure handle
      assert(isgraphics(figures{fig_i}{fig_j},'Figure'),'Invalid input type: choices are filename or figure handle');
      subplotFilename = sprintf('%s%ssubplot_%d_%d.png',dir,filesep,fig_i,fig_j);
      set(0, 'CurrentFigure', figures{fig_i}{fig_j});
      set(figures{fig_i}{fig_j}, 'color', pageColor);
      export_fig(subplotFilename,'-m1.2','-opengl');
      close(gcf);
    elseif ~isempty(regexp(figures{fig_i}{fig_j},'.fig','ONCE'))  %it's a .fig
      h = openfig(figures{fig_i}{fig_j},'invisible');
      set(gcf, 'color', pageColor);
      subplotFilename = regexprep(figures{fig_i}{fig_j},'.fig','.png');
      export_fig(subplotFilename,'-m1.2','-opengl');
      close(h);
    else %its a png
      subplotFilename = figures{fig_i}{fig_j};
    end
    
    %Load the .png
    pixels = imread(subplotFilename);
    figurePixels{fig_i}{fig_j} = pixels;
    [verticalSizes(fig_i,fig_j), horizontalSizes(fig_i,fig_j), ~] = size(pixels);
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
      if all(figures{fig_i+movesDown}{fig_j} == figures{fig_i}{fig_j}) %Line here has issue with [2;2;1] figure setup.
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
    delete(subplotFilename);
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
function [  ] = texSubplot( title, dir, fname, figures )
%texSubplot uses the beamer class of LaTeX to make a subplot figure
%   If same fig on sequential rows, fig spans those rows
%   Params: 
%   - title, string (possibly empty), title for the new figure
%   - dir, string, location at which to store tex file (omit trailing '/')
%   - fname: string, name for the tex file (omit the .tex extension)
%   - figures: cell array of cell arrays of figure handles, .fig filenames, or .png
%              filenames. The first cell array index gives the row index in the assembled figure; 
%               the second gives the column in the assembled final figure

% preprocess the figure inputs, if not already png files
for fig_i = 1:length(figures)
  for fig_j = 1:length(figures{1})
    % if a subplot is .fig or fig handle, generate png and update figures with its filename
    if ~ischar(figures{fig_i}{fig_j}) || ~isempty(regexp(figures{fig_i}{fig_j},'.fig','ONCE'))
      if ~ischar(figures{fig_i}{fig_j})  %it's a figure handle
        assert(isgraphics(figures{fig_i}{fig_j},'Figure'),'Invalid input type: choices are filename or figure handle');
        subplotFilename = sprintf('%s/subplot_%d_%d.png',dir,fig_i,fig_j);
        savefig(figures{fig_i}{fig_j},subplotFilename);
      else  %it's a .fig
        h = openfig(figures{fig_i}{fig_j});
        subplotFilename = regexprep(figures{fig_i}{fig_j},'.fig','.png');
        savefig(h,subplotFilename);
      end
      % now, replace all the entries for this figure with the filename
      movesRight = 1;
      while fig_j+movesRight < length(figures{1})
        if all(figures{fig_i}{fig_j+movesRight} == figures{fig_i}{fig_j}) %note: 'all' makes it work for strings and handles
          figures{fig_i}{fig_j+movesRight} = subplotFilename;
          movesRight = movesRight + 1;
        else
          break
        end
      end
      movesDown = 1;
      while fig_i+movesDown < length(figures)
        if all(figures{fig_i+movesDown}{fig_j} == figures{fig_i}{fig_j})
          for moveRight = 0:movesRight-1
            figures{fig_i+movesDown}{fig_j+moveRight} = subplotFilename;
          end
          movesDown = movesDown + 1;
        else
          break
        end
      end
      figures{fig_i}{fig_j} = subplotFilename;
    end
  end
end

mkdir(dir);
system(sprintf('cp subplotTemplate.txt %s',dir));
f = fopen(strcat(dir,'/','subplotTemplate.txt'),'a');
fprintf(f,'\\begin{frame}\n');
if ~isempty(title)
  fprintf(f,'\\frametitle{%s}\n', title);
end
fprintf(f,'\\begin{columns}[t]\n');
for column_i = 1:length(figures{1})
  fprintf(f,'\\column{%.2f\\textwidth}\n',1/size(figures,2));
%   rowSpan = 1;
  for row_i = 1:length(figures)
%     if row_i < length(figures) && strcmp(figures{row_i}{column_i}, figures{row_i+1}{column_i})
%       rowSpan = rowSpan + 1;
%     else
%       fprintf(f,'\\includegraphics[width=%.2f\\textwidth,keepaspectratio]{%s}\n',.8/length(figures{1}),figures{row_i}{column_i});
%       %fprintf(f,'\\includegraphics[height=%.2f\\textheight,keepaspectratio]{%s}\n',rowSpan/length(figures),figures{row_i}{column_i});
%       rowSpan = 1;
%     end
    
    fprintf(f,'\\includegraphics[width=%.2fcm,height=%.2fcm,keepaspectratio]{%s}\\\\\n',8/length(figures{1}),6.5/length(figures),figures{row_i}{column_i});
    fprintf(f,'\\vspace{0.3cm}\n');
  end
end
fprintf(f,'\\end{columns}\n\\end{frame}\n\\end{document}');
fclose(f);
system(sprintf('mv %s/subplotTemplate.txt %s/%s.tex',dir,dir,fname));
if isempty(regexp(getenv('PATH'),'TeX/texbin','ONCE'))
  setenv('PATH',strcat(getenv('PATH'),':/usr/texbin:/Library/TeX/texbin'));
end
wkdir = pwd();
cd(dir);
system(sprintf('pdflatex %s/%s',dir,fname));
cd(wkdir);
end


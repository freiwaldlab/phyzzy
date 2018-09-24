function [  ] = texSubplot( title, dir, fname, figures )
%texSubplot uses the beamer class of LaTeX to make a subplot figure
%   If same fig on sequential rows, fig spans those rows
%   Params: 
%   - title, string (possibly empty), title for the new figure
%   - dir, string, location at which to store tex file (omit trailing '/')
%   - fname: string, name for the tex file (omit the .tex extension)
%   - figures: cell array of cell arrays of figure handles, .fig filenames, or .png
%              filenames (currently, only .png supported). The first cell
%              array index gives the row index in the assembled figure; the second
%              gives the column in the assembled final figure

% preprocess the figure inputs, if not already png files
for fig_i = 1:length(figures)
  for fig_j = 1:length(figures{1})
    if ~ischar(figures{fig_i}{fig_j}) || isempty(regexp(figures{fig_i}{fig_j},'.png','ONCE'))
      %todo: implement...
      disp('Support for .fig and figure handles not yet implemented');
      return;
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


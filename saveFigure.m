function [ ] = saveFigure( outDir, filename, figData, saveFig, exportFig, saveData, varargin )
%UNTITLED11 Summary of this function goes here
%   Detailed explanation goes here

if ~isempty(varargin)
  footer = varargin{1};
  ax = axes('Position',[0 0 1 .05], 'Visible','off');
  text(ax, .025,.5,footer, 'fontsize',12);
end
if saveFig     
  savefig(strcat(outDir, filename,'.fig'));
end
if exportFig
  export_fig([outDir strcat(filename, '.png')],'-m1.2','-transparent','-opengl');
end
if saveData
  save(strcat(outDir,filename,'_data.mat'),'figData');
end
if length(varargin) > 1 && varargin{2} == 'close'
  close()
end
end


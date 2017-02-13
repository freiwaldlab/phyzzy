function draw_vert_line(x,varargin)
% Draw a vertical line, for the whole visible
% y axis range
syl = ylim;
if exist('color','var')
    line([x x],[syl(1) syl(2)],varargin{:});
else
    line([x x],[syl(1) syl(2)],varargin{:});
end

function lineHand = draw_vert_line(x,varargin)
% Draw a vertical line, for the whole visible.
% y axis range

if exist('color','var')
    lineHand = line([x x], ylim(),varargin{:});
else
    lineHand = line([x x], ylim(),varargin{:});
end

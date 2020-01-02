function OutCellArray = initNestedCellArray(varargin)
% function initializes an empty nested cell array of a particular tiered
% structure. First Input can either be
% - an array, where each number represents the max index at that tier
% OR
% - a nested cell array, which will be copied in structure
% 2nd Argument = 'ones','zeros','NaN', or 'cell', to distinguish structure
% at bottom. default is cell
% 3rd argument = array of dimensions to feed into 2nd argument function.
% default is [1, 1].
% 4th argument = depth, if the full depth isn't desired. (default = 100)
%% Parse Arguments
fill = 'cell';
fillSize = [1 1];
depth = 100;

switch nargin()
  case 0
    error('Not enough inputs, minimum 1')
  case 1
    inCell = varargin{:};
  case 2
    [inCell, fill] = varargin{:};
  case 3
    [inCell, fill, fillSize] = varargin{:};
  case 4
    [inCell, fill, fillSize, depth] = varargin{:};
  case 5  % meant to allow for recursive functioning of cells with doubles.
    [inCell, fill, fillSize, depth, originalInput] = varargin{:};
end

if nargin ~= 5
  originalInput = class(inCell);
end


%% Function
%To-do - use this function to initialize all matching size and format
%arrays above.

if depth == 0 || ~strcmp(class(inCell), originalInput) || isempty(inCell)
  switch fill
    case 'zeros'
      OutCellArray = zeros(fillSize);
    case 'ones'
      OutCellArray = ones(fillSize);
    case 'NaN'
      OutCellArray = NaN(fillSize);
    case 'cell'
      OutCellArray = cell(fillSize);
  end
else
  switch originalInput
    case 'cell'
      OutCellArray = cell(size(inCell));
      for inCell_ind = 1:length(inCell)
        OutCellArray{inCell_ind} = initNestedCellArray(inCell{inCell_ind}, fill, fillSize, depth-1, originalInput);
      end
    case 'double'
      OutCellArray = cell(1,inCell(1));
      for inCell_ind = 1:length(OutCellArray)
        OutCellArray{inCell_ind} = initNestedCellArray(inCell(2:end), fill, fillSize, depth-1, originalInput);
      end
  end
end
end

function OutCellArray = fillOutput(fill, fillSize)

end

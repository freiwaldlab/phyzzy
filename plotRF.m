function [ fh ] = plotRF( rfGrid, data, errors, method, colors )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

fh = figure();
switch method
  case 'stem'
    assert(size(data,1) == size(rfGrid,1) && size(data,2) == 1);
    stem3(rfGrid(:,1), rfGrid(:,2),data,'linewidth',10);
    hold on
    if ~isempty(errors)
      for grid_i = 1:size(rfGrid,1)
        plot3([rfGrid(grid_i,1), rfGrid(grid_i,1)], [rfGrid(grid_i,2), rfGrid(grid_i,2)],[data(grid_i) - errors(grid_i), data(grid_i) + errors(grid_i)],'marker','none','color','k');
      end
    end
  case 'bar'
    disp('RF bar plot not currently implemented');
    %fcn to use is bar3
    return
  case 'subplots'
    if length(unique(rfGrid(1,:))) == 2*length(rfGrid(rfGrid(1,:) == rfGrid(1,1))) % triangular grid, todo: also check for y vals to get triangle vs hexagon
      numX = length(unique(rfGrid(1,:))); %note BUG!!!: indexing is backwards here
      numY=  length(unique(rfGrid(2,:)));
      gridXsorted = sort(rfGrid(1,:));
      gridYsorted = sort(rfGrid(2,:));
      for grid_i = 1:size(rfGrid,1)
        pt = rfGrid(grid_i); %pt stands for 'point'
        xInd = find(gridXsorted == pt(1));
        yInd = find(gridYsorted == pt(2));
        subplots(numX,numY,xInd+(yInd-1)*numX:xInd+1+(yInd-1)*numX);
        if size(data(grid_i),2) == 1
          superbar(data(grid_i),'E',errors(grid_i));
        else
          mseb(data(grid_i,:),errors(grid_i,:));
        end
      end
    elseif length(unique(rfGrid(1,:))) == length(rfGrid(rfGrid(1,:) == rfGrid(1,1))) %square grid
    else
      disp('subplots RF plot implemented only for square and triangular grids');
      return
    end
  case 'heat'
    numX = length(unique(rfGrid(:,1)));
    numY=  length(unique(rfGrid(:,2)));
    gridXsorted = unique(sort(rfGrid(:,1)));
    gridYsorted = unique(sort(rfGrid(:,2)));
    dataMap = NaN(numX+1,2*numY);
    errorMap = NaN(numX+1,2*numY);
    % first, build the dataMap with (x,y) indexing; we will need to
    % transpose it to display it correctly as an image
    for grid_i = 1:size(rfGrid,1)
      pt = rfGrid(grid_i,:); %pt stands for 'point'
      xInd = find(gridXsorted == pt(1));
      yInd = find(gridYsorted == pt(2));
      dataMap(xInd:xInd+1,2*yInd-1:2*yInd) = data(grid_i);
      errorMap(xInd:xInd+1,2*yInd-1:2*yInd) = errors(grid_i);
    end
    subplot(2,3,[1,2,4,5]);
    imAlpha = ~isnan(dataMap');
    imagesc([min(rfGrid(:,2)) max(rfGrid(:,2))], [min(rfGrid(:,1)) max(rfGrid(:,1))], dataMap','AlphaData',imAlpha);
    set(gca,'YDir','normal');
    set(gca,'color',.5*[1 1 1]); %makes NaN positions grey
    set(gca,'DataAspectRatio',[1 1 1]);
    ylabel('vertical position (dva)');
    xlabel('horizontal position (dva)');
    colormap(colors);
    caxis([min(data-errors) max(data+errors)]);
    h = colorbar();
    ylabel(h,'Hz','FontSize',14);
    
    subplot(2,3,3);
    imagesc([min(rfGrid(:,2)) max(rfGrid(:,2))], [min(rfGrid(:,1)) max(rfGrid(:,1))],dataMap'+errorMap','AlphaData',imAlpha);
    set(gca,'YDir','normal');
    set(gca,'color',.5*[1 1 1]); %makes NaN positions grey
    set(gca,'DataAspectRatio',[1 1 1]);
    ylabel('vertical position (dva)');
    xlabel('horizontal position (dva)');
    colormap(colors);
    caxis([min(data-errors) max(data+errors)]);
    
    subplot(2,3,6);
    imagesc([min(rfGrid(:,2)) max(rfGrid(:,2))], [min(rfGrid(:,1)) max(rfGrid(:,1))],dataMap'-errorMap','AlphaData',imAlpha);
    set(gca,'YDir','normal');
    set(gca,'color',.5*[1 1 1]); %makes NaN positions grey
    set(gca,'DataAspectRatio',[1 1 1]);
    ylabel('vertical position (dva)');
    xlabel('horizontal position (dva)');
    colormap(colors);
    caxis([min(data-errors) max(data+errors)]);
    
  case 'cnr'  %plot (value - mean) / sem; todo: cnr relative to baseline
    numX = length(unique(rfGrid(:,1)));
    numY=  length(unique(rfGrid(:,2)));
    gridXsorted = unique(sort(rfGrid(:,1)));
    gridYsorted = unique(sort(rfGrid(:,2)));
    dataMap = NaN(numX+1,2*numY);
    errorMap = NaN(numX+1,2*numY);
    % first, build the dataMap with (x,y) indexing; we will need to
    % transpose it to display it correctly as an image
    for grid_i = 1:size(rfGrid,1)
      pt = rfGrid(grid_i,:); %pt stands for 'point'
      xInd = find(gridXsorted == pt(1));
      yInd = find(gridYsorted == pt(2));
      dataMap(xInd:xInd+1,2*yInd-1:2*yInd) = data(grid_i);
      errorMap(xInd:xInd+1,2*yInd-1:2*yInd) = errors(grid_i);
    end
    cnrMap = (dataMap - mean(mean(dataMap(~isnan(dataMap))))) ./ errorMap;
    imAlpha = ~isnan(cnrMap');
    imagesc([min(rfGrid(:,2)) max(rfGrid(:,2))], [min(rfGrid(:,1)) max(rfGrid(:,1))], cnrMap','AlphaData',imAlpha);
    set(gca,'YDir','normal');
    set(gca,'color',.5*[1 1 1]);  %makes NaN positions grey
    set(gca,'DataAspectRatio',[1 1 1]);
    ylabel('vertical position (dva)');
    xlabel('horizontal position (dva)');
    %colormap(colors);
    h = colorbar();
    ylabel(h,'CNR','FontSize',14);
     
  case 'modulation' %the idea here is to show the significance of the NN and 2ndNN differences from uniform, centered at each point
    
  case 'surf'
    assert(size(data,1) == size(rfGrid,1) && size(data,2) == 1);
    [X,Y] = meshgrid(sort(unique(rfGrid(:,1))),sort(unique(rfGrid(:,2))));
    dataGrid = griddata(rfGrid(:,1),rfGrid(:,2),data,X,Y,'cubic');
    surf(X,Y,dataGrid,'LineStyle','-','FaceAlpha','0.25','EdgeColor','interp','FaceColor','interp');
    colormap(colors);
    xlabel('horizontal position (dva)');
    ylabel('vertical position (dva)');
    zlabel('firing rate (Hz)');
    axis tight
    colorbar();
    hold on
    plot3(rfGrid(:,1), rfGrid(:,2),data,'marker','o','color','k','linestyle','none');
    if ~isempty(errors)
      for grid_i = 1:size(rfGrid,1)
        plot3([rfGrid(grid_i,1), rfGrid(grid_i,1)], [rfGrid(grid_i,2), rfGrid(grid_i,2)],[data(grid_i), data(grid_i) + errors(grid_i)],'marker','none','color','k');
      end
    end
  case 'splitHalf'
    disp('split half rf plot not yet implemented');
  case 'fisherInfo'
    disp('fisher info rf plot not yet implemented');
end
end


function display_map(x,y,z,xi,yi,sgaus,normalize,save,figTitle,filename)
%x,y is position in 2D map, z is assigned value for each point. xi,yi is
%interpolated x/y axis. sgaus (sigma of gaussian envelope at each point)
%is a parameter for interpolation which empirically determined. Michael
%used 2.2857*gridsize. normalize set to zero
%sgaus = 2.2857;

lambdas = RBFcoef(z, [x y]', sgaus);

[mxi,myi] = meshgrid(xi,yi);
zi = reshape( ...
	RBFipol(lambdas, [x y]', [mxi(:) myi(:)]', sgaus, normalize), ...
	size(mxi));
figure();
imagesc(xi,yi,zi);
set(gca,'YDir','normal');
set(gca,'DataAspectRatio',[1 1 1]);
ylabel('vertical position (dva)');
xlabel('horizontal position (dva)');
hold on;
plot(x,y,'wo');
hold off;
load AkinoriColormap.mat; colormap(AkinoriColormap);
colorbar;
title(figTitle);
if save
  export_fig(filename,'-m1.2','-transparent','-opengl'); 
end
end



%% Convert a NEV file into a MAT file for fast data read and smaller file
%  size
%
%  Use NEV2MAT
%
%  Version 1.1
%  September 02, 2009

function NEV2MAT

datafolder = uigetdir('D:\', 'Select the folder containing Cerebus folders...');
savefolder = uigetdir('D:\', 'Select the folder where you would like to save MAT files...');
foldernames = ls(datafolder);
foldercounts = size(foldernames,1);
cd(datafolder)

for i=3:foldercounts-1
    clear NEV;
    foldername = foldernames(i,1:end);
    filenames = ls(foldername);
    trialCount = max(str2num(filenames(1:end-4,end-4)));
    for ii = 1:trialCount
        fullpath = [datafolder '\' foldername '\' foldername '-00' num2str(ii) '.nev'];
        NEV = openNEV3('report', fullpath);
    
        nfullpath = [savefolder '\' foldername '-00' num2str(ii) '.mat'];

        save(nfullpath, 'NEV');
    end
end
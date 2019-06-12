function behaviorsummaryPhyzzy(filename)
% NIMH MonkeyLogic
%
% This function reads NIMH ML data files (*.bhv2; *.h5; *.mat) and shows a
% performance summary. Version modified to work with Phyzzy, and present
% additional data.

% select file
if exist('filename','var')
    if 2~=exist(filename,'file'), return, end
else
  [file,path] = uigetfile({'*.mat;*.bhv2;*.h5;','MonkeyLogic Files (*.mat,*.bhv2,*.h5)'},...
   'Select a File');
  filename = [path file];
end

%Load the necessary data and set the appropriate paths
[data,MLConfig,TrialRecord,datafile] = mlread(filename);
[filepath,n,e] = fileparts(datafile);
filename = [n e];

% Unpack data
trial = [data.Trial];
block = [data.Block];
blockswitch = find(diff(block));
blockorder = block([1 blockswitch+1]);
errortype = [data.TrialError];
ntrial = length(trial);
corrInd = find(errortype == 0, 1);
userVars = data(corrInd).UserVars; %Created by bhv_variable fxn.
userVars2 = data(corrInd).VariableChanges; % Created by editable fxn.

%
filebits = strsplit(filepath, '/');
filebits = cell2mat(fullfile(filebits(1:end-1),'/'));
%If called individually on PC
if isempty(filebits)
  filebits = strsplit(filepath, '\');
  filebits = cell2mat(fullfile(filebits(1:end-1),'\'));
end
xlsxLog = dir([filebits 'Recordings*.xlsx']);
%If you can find an xls log file, open it.
if ~isempty(xlsxLog)
  %Find the index for data based on the titles of the columns
  [~, ~, raw] = xlsread([xlsxLog.folder filesep xlsxLog.name]);
  recID = strcmp(raw(1,:),'recording');
  rowID = find(strcmp(raw(:,recID),n));
  %If there is a matching recordingID, start pulling variables from it.
  if ~isempty(rowID) %If the row is present, use any variables present to overwrite what was present in the MKL Table +
    %(for instances b/t July - Dec where runs weren't updated, correct info is pulled from Log book).
    %Pull the data from those columns, using the "Recording" column to
    %match days to runs.
    electrodeInd = strcmp(raw(1,:),'impedance (MOhms)');
    holeMLInd = find(strcmp(raw(1,:),'ML'));
    holeAPInd = find(strcmp(raw(1,:),'AP'));
    recordingInd = strcmp(raw(1,:),'putative distance');
    recActComInd = strcmp(raw(1,:),'comments on signal');
    regionInd = strcmp(raw(1,:),'putative region');
    commInd = strcmp(raw(1,:),'online analysis');
    
    userVars2Temp.activityComments = raw{rowID,recActComInd};
    userVars2Temp.comments = raw{rowID,commInd};
    userVars2Temp.electrodeImp = raw{rowID,electrodeInd};
    userVars2Temp.gridHole = [num2str(raw{rowID,holeMLInd}) '-' num2str(raw{rowID,holeAPInd})];
    userVars2Temp.putativeRegion = raw{rowID,regionInd};
    userVars2Temp.recordingDepth = raw{rowID,recordingInd};
    userVars2Temp.rewardFreq = 1;
    
    %Extract values from XLS are sometimes NaN, which must be treated
    %correctly.
    varFields = fields(userVars2Temp);
    for field_i = 1:length(fields(userVars2Temp))
      %If the XLS has anything, we assume its right, if the experimental
      %data struct is missing the field, or has it empty, we put NaN
      %anyway.
      if any(~isnan(eval(sprintf('userVars2Temp.%s', varFields{field_i})))) || (isempty(eval(sprintf('userVars2.%s', varFields{field_i}))) || ~isfield(userVars2, varFields{field_i}))
          eval(sprintf('userVars2.%s = userVars2Temp.%s;', varFields{field_i}, varFields{field_i}))
      end
    end
  end
end

%If we're missing other reward related variables, which weren't logged in
%Aug/July, fill them in from the trial data.
if ~isfield(userVars,'clipDuration')
  userVars.clipStartTime = 0;
  userVars.clipDuration = data(corrInd).ObjectStatusRecord.SceneParam(2).AdapterArgs{1,3}{3};
  RewardTmp = [data.RewardRecord];
  userVars.rewardTimeMean = round(mean([RewardTmp.EndTimes]-[RewardTmp.StartTimes]));
  userVars.rewardSD = 0; %Experiments where this isn't recorded were before this was used, and therefore must be close to 0.
end

% constants
colororder = [0 1 0; 0 1 1; 1 1 0; 0 0 1; 0.5 0.5 0.5; 1 0 1; 1 0 0; .3 .7 .5; .7 .2 .5; .5 .5 1; .75 .75 .5];
corder(1,1:11,1:3) = colororder(1:11,1:3);
%figure_bgcolor = [.65 .70 .80];

% create figure
% fw = 800;
% fh = 600;
% screen_pos = GetMonitorPosition(mglgetcommandwindowrect);
% h = findobj('tag','mlmainmenu');
% if isempty(h), pos = screen_pos; else pos = get(h,'position'); end
% fx = pos(1) + 0.5 * (pos(3) - fw);
% if fx < screen_pos(1), fx = screen_pos(1) + 8; end
% fy = min(pos(2) + 0.5 * (pos(4) - fh),screen_pos(2) + screen_pos(4) - fh - 85);
% fig_pos = [fx fy fw fh];
hFig =   figure('Name','MonkeyLogic Behavioral Summary','NumberTitle','off');
set(hFig, 'tag','behaviorsummary', 'numbertitle','off'); %'color',figure_bgcolor);

% performance plot
perfplot = subplot('position',[0 0.5 1 0.5]);
set(gca, 'ylim',[0 1], 'position',[0.08 0.65 0.5 0.3], 'box','on');

hold on;
h = title(filename,'interpreter','none');
set(h,'fontsize',10,'fontweight','bold');
xlabel('Trial number');
ylabel('Fraction correct');

smoothwin = 10;
if ntrial < 5*smoothwin  % histogram when the number of trials is <50
  x = 1:ntrial;
  for m=0:9
    y = double(m==errortype);
    h = bar(x,y,1);
    set(h,'facecolor',colororder(m+1,:),'edgecolor',[1 1 1]);
  end
  set(gca,'xlim',[0.5 ntrial+0.5],'xtick',1:ntrial);
else                     % smoothe the performance curve if there are 50 trials or more
  kernel = [0.0047 0.0087 0.0151 0.0245 0.0371 0.0525 0.0693 0.0853 0.0979 0.1050 ...
    0.1050 0.0979 0.0853 0.0693 0.0525 0.0371 0.0245 0.0151 0.0087 0.0047]; % gaussian window
  yarray1 = zeros(ntrial, 12);
  for m = 0:10
    r = conv(double(m==errortype),kernel,'same')';
    yarray1(:,m+2) = yarray1(:,m+1) + r;
  end
  xarray1 = (1:ntrial)';
  xarray1 = repmat(xarray1,1,12);
  xarray2 = flipud(xarray1);
  yarray2 = flipud(yarray1);
  x = cat(1,xarray1(:,1:11),xarray2(:,2:12));
  y = cat(1,yarray1(:,1:11),yarray2(:,2:12));
  patch(x, y, corder);
  set(gca,'xlim',[1 ntrial],'ytick',0:0.25:1);
  
  hline(1) = line([0 ntrial],[0.5 0.5]);
  hline(2) = line([0 ntrial],[0.25 0.25]);
  hline(3) = line([0 ntrial],[0.75 0.75]);
  set(hline,'color',[0.7 0.7 0.7],'linewidth',1);
  
  nblock = length(blockswitch);
  h = zeros(nblock,1);
  ht = h;
  texty = 0.05;
  for m = 1:nblock
    x1 = blockswitch(m);
    h(m) = line([x1 x1], [0 1]);
    if 1<m, x2 = blockswitch(m-1); else x2 = 0; end
    xm = (x1 + x2)/2;
    ht(m) = text(xm, texty, num2str(blockorder(m)));
  end
  if ~isempty(h)
    xm = (blockswitch(m) + length(trial))/2;
    ht(m+1) = text(xm, texty, num2str(blockorder(m+1)));
    set(h,'color',[1 1 1],'linewidth',1);
  else
    xm = ntrial/2;
    ht = text(xm, texty, num2str(blockorder));
  end
  set(ht, 'horizontalalignment','center', 'color',[1 1 1], 'fontweight','bold', 'fontsize',14);
end

% Recording Site Image
recordingSite = subplot('position',[0.6 0.65 0.35 0.3]);
monthYearRec = data(1).TrialDateTime(1)*100+data(1).TrialDateTime(2); %Determine the month/year of the recording
if monthYearRec < 201808
  %If recording took place prior to August 2018, Old Grid
  gridPic = imread('C:\Data 2018\GridPics\Grid_v1.png');
  gridOrigin = [37 159];
  gridStep = 34;
  holeDiam = 14;
  MLOffset = 1;
else %monthYearRec < 201901 %New Grid was used (August - Dec 2018)
  gridPic = imread('C:\Data 2018\GridPics\Grid_v2.png');
  gridOrigin = [60 70];
  gridStep = 45.5;
  holeDiam = 19;
  MLOffset = 4;
end

% Look for the number of holes and plot each correctly.
holes = split(userVars2.gridHole,[";",","]);
if length(holes) > 1
  removeSpaces = @(cellX) cellX(~isspace(cellX));
  holes = cellfun(removeSpaces, holes,'UniformOutput',false);
end

%Draw each of the holes, appropriately given the way they are described
for hole_i = 1:length(holes)
  if isletter(userVars2.gridHole(1))
    APGridHole = str2double(holes{hole_i}(2:end))-1;
    letters = {'A','B','C','D','E'};
    MLGridHole = find(strcmp(letters, holes{hole_i}(1)))+3;
    holeCords = [APGridHole (MLGridHole-MLOffset)];
  else
    APGridHole = (str2double(holes{hole_i}(3:end))-1);
    holeCords = [APGridHole (str2double(holes{hole_i}(1))-MLOffset)];
  end
  recHoleSite = [gridOrigin + holeCords.*gridStep holeDiam];
  gridPic = insertShape(gridPic,'FilledCircle',recHoleSite,'color','black','Opacity',1);
end

imagesc(gridPic);
recordingSite.YTickLabel = [];
recordingSite.XTickLabel = [];
title('Grid and Recording Sites')

% task info
x = 45; y = 220; w = 270; lineinterval = 19; fontsize = 10;
uicontrol('parent',hFig,'style','text','position',[x y w-60 25],'string',['Exp: ' MLConfig.ExperimentName], 'fontsize',fontsize,'fontweight','bold','horizontalalignment','left');

%,'units','pixel'
if isfield(MLConfig,'MLPath') && isfield(MLConfig.MLPath, 'ConditionsFile')
  [~,Conditions] = fileparts(MLConfig.MLPath.ConditionsFile);
  y = y - lineinterval;
  uicontrol('parent',hFig,'style','text','units','pixel','position',[x y w 22],'string',['Condition: ' Conditions], 'fontsize',fontsize,'fontweight','bold','horizontalalignment','left');
end

y = y - lineinterval;
uicontrol('parent',hFig,'style','text','units','pixel','position',[x y w 22],'string',['Subject: ' MLConfig.SubjectName], 'fontsize',fontsize,'fontweight','bold','horizontalalignment','left');

if ~isempty(TrialRecord) && isfield(TrialRecord,'TaskInfo')
  StartTime = TrialRecord.TaskInfo.StartTime;
  EndTime = TrialRecord.TaskInfo.EndTime;
else
  StartTime = datenum(data(1).TrialDateTime);
  EndTime = datenum(data(end).TrialDateTime) + data(end).BehavioralCodes.CodeTimes(end)/86400000;
end

% Write information below figures
y = y - lineinterval;
uicontrol('parent',hFig,'style','text','units','pixel','position',[x y w 22],'string',['Start time: ' datestr(StartTime)], 'fontsize',fontsize,'fontweight','bold','horizontalalignment','left');

y = y - lineinterval;
uicontrol('parent',hFig,'style','text','units','pixel','position',[x y w 22],'string',['End time: ' datestr(EndTime)], 'fontsize',fontsize,'fontweight','bold','horizontalalignment','left');

y = y - lineinterval;
uicontrol('parent',hFig,'style','text','units','pixel','position',[x+80 y w-80 22],'string',sprintf('(Elapsed: %s)',datestr(EndTime-StartTime,'HH:MM:SS')), 'fontsize',fontsize,'fontweight','bold','horizontalalignment','left');

y = y - lineinterval;
uicontrol('parent',hFig,'style','text','units','pixel','position',[x y w 22],'string',sprintf('Stimulus duration: %d ms',userVars.clipDuration), 'fontsize',fontsize,'fontweight','bold','horizontalalignment','left');

y = y - lineinterval;
uicontrol('parent',hFig,'style','text','units','pixel','position',[x y w 22],'string',sprintf('Reward amount (mean): %d ms',userVars.rewardTimeMean), 'fontsize',fontsize,'fontweight','bold','horizontalalignment','left');

y = y - lineinterval;
uicontrol('parent',hFig,'style','text','units','pixel','position',[x y w 22],'string',sprintf('Reward SD: %d ms',userVars.rewardSD), 'fontsize',fontsize,'fontweight','bold','horizontalalignment','left');

y = y - lineinterval;
nblock = length(blockorder);
uicontrol('parent',hFig,'style','text','units','pixel','position',[x y w 22],'string',sprintf('# of completed blocks: %d',nblock), 'fontsize',fontsize,'fontweight','bold','horizontalalignment','left');

correct = sum(0==errortype);
y = y - lineinterval;
uicontrol('parent',hFig,'style','text','units','pixel','position',[x y w 22],'string',sprintf('Success rate: %.2f%% (= %d / %d)',correct/ntrial*100,correct,ntrial), 'fontsize',fontsize,'fontweight','bold','horizontalalignment','left');

attempted = sum(errortype ~= 4);
y = y - lineinterval;
uicontrol('parent',hFig,'style','text','units','pixel','position',[x y w 22],'string',sprintf('True success rate: %.2f%% (= %d / %d)',correct/attempted*100,correct,attempted), 'fontsize',fontsize,'fontweight','bold','horizontalalignment','left');

%Draw to the next column;
lineinterval = 21;

x = 305; y = 210; w = 270;
uicontrol('parent',hFig,'style','text','units','pixel','position',[x y w 22],'string',sprintf('Grid Hole: %s   Region: %s',userVars2.gridHole, userVars2.putativeRegion), 'fontsize',fontsize,'fontweight','bold','horizontalalignment','left');

y = y - lineinterval;
uicontrol('parent',hFig,'style','text','units','pixel','position',[x y w 22],'string',sprintf('Electrode (MOhms): %s',num2str(userVars2.electrodeImp)), 'fontsize',fontsize,'fontweight','bold','horizontalalignment','left');

y = y - lineinterval;
uicontrol('parent',hFig,'style','text','units','pixel','position',[x y w 22],'string',sprintf('Recording Depth: %s um',num2str(userVars2.recordingDepth)), 'fontsize',fontsize,'fontweight','bold','horizontalalignment','left');

y = y - lineinterval;
uicontrol('parent',hFig,'style','text','units','pixel','position',[x y w 22],'string',sprintf('Reward Freq: %s',num2str(userVars2.rewardFreq)), 'fontsize',fontsize,'fontweight','bold','horizontalalignment','left');

y = y - lineinterval*2.5;
uicontrol('parent',hFig,'style','text','units','pixel','position',[x y w-20 22*2.5],'string',sprintf('Activity comments: %s',userVars2.activityComments), 'fontsize',fontsize,'fontweight','bold','horizontalalignment','left');

if isfield(userVars2,'comments')
  y = y - lineinterval*3;
  uicontrol('parent',hFig,'style','text','units','pixel','position',[x y w-20 22*3],'string',sprintf('General Comments: %s',userVars2.comments), 'fontsize',fontsize,'fontweight','bold','horizontalalignment','left');
end

end

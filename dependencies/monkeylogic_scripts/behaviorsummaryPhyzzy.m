function behaviorsummaryPhyzzy(filename)
% NIMH MonkeyLogic
%
% This function reads NIMH ML data files (*.bhv2; *.h5; *.mat) and shows a
% performance summary.

%Stripped down to present less info and make compatible with phyzzy.

% select file
if exist('filename','var')
    if 2~=exist(filename,'file'), return, end
else
    filename = '';
end
[data,MLConfig,TrialRecord,datafile] = mlread(filename);
[filepath,n,e] = fileparts(datafile);
filename = [n e];

% data
trial = [data.Trial];
block = [data.Block];
blockswitch = find(diff(block));
blockorder = block([1 blockswitch+1]);
errortype = [data.TrialError];
ntrial = length(trial);
userVars = data(end).UserVars;
userVars2 = data(end).VariableChanges;

if ~isfield(userVars,'clipDuration')
    filebits = strsplit(filepath, '/');
    filebits = cell2mat(fullfile(filebits(1:end-1),'/'));
    xlsxLog = dir([filebits 'Recordings*.xlsx']);
    [num, txt, raw] = xlsread([xlsxLog.folder filesep xlsxLog.name]);
    
    %Find the index for data based on the titles of the columns
    electrodeInd = find(strcmp(raw(1,:),'impedance (MOhms)'));
    holeMLInd = find(strcmp(raw(1,:),'ML'));
    holeAPInd = find(strcmp(raw(1,:),'AP'));
    recordingInd = find(strcmp(raw(1,:),'putative distance'));
    thresInd = find(strcmp(raw(1,:),'threshold (microVolts)'));
    recActComInd = find(strcmp(raw(1,:),'comments on signal'));
    recID = find(strcmp(raw(1,:),'recording'));
    regionInd = find(strcmp(raw(1,:),'putative region'));
    commInd = find(strcmp(raw(1,:),'online analysis'));
    rowID = find(strcmp(raw(:,recID),n));
    
    %Pull the data from those columns, using the "Recording" column to
    %match days to runs.
    try
        userVars.clipDuration = num2str(data(1).ObjectStatusRecord.SceneParam(2).AdapterArgs{1,3}{3});
    catch
        userVars.clipDuration = num2str(data(10).ObjectStatusRecord.SceneParam(2).AdapterArgs{3}{1,2});
    end
    userVars.rewardSD = '0';
    userVars.rewardTimeMean = MLConfig.RewardFuncArgs.Duration;
    userVars2.comments = num2str(raw{rowID,commInd});
    userVars2.gridHole = [ num2str(raw{rowID,holeMLInd}) '-' num2str(raw{rowID,holeAPInd})]; %#ok<FNDSB>
    userVars2.electrodeImp = num2str(raw{rowID,electrodeInd});
    userVars2.recordingDepth = num2str(raw{rowID,recordingInd});
    userVars2.putativeRegion = num2str(raw{rowID,regionInd});
    userVars2.thresholdSetting = num2str(raw{rowID,thresInd});
    userVars2.activityComments = num2str(raw{rowID,recActComInd});
    
end

% constants
colororder = [0 1 0; 0 1 1; 1 1 0; 0 0 1; 0.5 0.5 0.5; 1 0 1; 1 0 0; .3 .7 .5; .7 .2 .5; .5 .5 1; .75 .75 .5];
corder(1,1:11,1:3) = colororder(1:11,1:3);
figure_bgcolor = [.65 .70 .80];

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
hFig = figure;
set(hFig, 'tag','behaviorsummary', 'numbertitle','off', 'name',filename); %'color',figure_bgcolor);

% performance plot
subplot('position',[0 0.5 1 0.5]);
set(gca, 'ylim',[0 1], 'position',[0.08 0.65 0.88 0.3], 'box','on');

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

% rt plot
% subplot('position',[0.5 0 0.5 0.5]);
% set(gca,'Units','normalized','position',[0.52 0.11 0.44 0.35], 'box','on', 'color',figure_bgcolor);
% 
% userplotfunc = [];
% if ~isempty(MLConfig.UserPlotFunction) && ~isempty(TrialRecord)
%     [~,n,e] = fileparts(MLConfig.UserPlotFunction);
%     userplotfunc = get_function_handle(MLConfig.UserPlotFunction);
%     if isempty(userplotfunc), userplotfunc = get_function_handle([filepath filesep n e]); end
%     if isempty(userplotfunc), userplotfunc = get_function_handle([pwd filesep n e]); end
%     if isempty(userplotfunc) && isfield(TrialRecord,'TaskInfo') && ~isempty(TrialRecord.TaskInfo.UserPlotFunction)
%         fid = fopen([tempdir n e],'w');
%         fwrite(fid,TrialRecord.TaskInfo.UserPlotFunction);
%         fclose(fid);
%         userplotfunc = get_function_handle([tempdir n e]);
%     end
%     if isempty(userplotfunc) && isfield(TrialRecord.User,'UserPlotFunction') && ~isempty(TrialRecord.User.UserPlotFunction)
%         fid = fopen([tempdir n e],'w');
%         fwrite(fid,TrialRecord.User.UserPlotFunction);
%         fclose(fid);
%         userplotfunc = get_function_handle([tempdir n e]);
%     end
%     try
%         if ~isempty(userplotfunc), userplotfunc(TrialRecord); end
%     catch
%         userplotfunc = [];
%     end
% end
% if isempty(userplotfunc)
%     reaction_time = [data.ReactionTime];
%     if any(~isnan(reaction_time))
%         [n,x] = hist(reaction_time,20);
%         h = bar(x,n,1);
%         set(h, 'facecolor',[1 1 0], 'edgecolor',[.8 .8 .8], 'linewidth',0.5);
%         set(gca, 'xlim',[0 max(x)*21/20], 'ylim',[0 1.1*max(n)]);
%     else
%         h = text(mean(get(gca,'xlim')), mean(get(gca,'ylim')), 'No Reaction Time Data Available');
%         set(h, 'fontsize',14, 'color',[.3 0 0], 'fontweight','bold', 'horizontalalignment','center');
%     end
%     xlabel('Reaction Time (msec)');
%     ylabel('Number of Trials');
% end

% task info
x = 45; y = 220; w = 270; lineinterval = 19; fontsize = 10;
uicontrol('parent',hFig,'style','text','position',[x y w-60 25],'string',['Experiment: ' MLConfig.ExperimentName], 'fontsize',fontsize,'fontweight','bold','horizontalalignment','left');

%,'units','pixel'
if isfield(MLConfig,'MLPath') && isfield(MLConfig.MLPath, 'ConditionsFile')
    [~,Conditions] = fileparts(MLConfig.MLPath.ConditionsFile);
    y = y - lineinterval;
    uicontrol('parent',hFig,'style','text','units','pixel','position',[x y w 22],'string',['Conditions: ' Conditions], 'fontsize',fontsize,'fontweight','bold','horizontalalignment','left');
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
uicontrol('parent',hFig,'style','text','units','pixel','position',[x y w 22],'string',sprintf('Electrode(MOhms): %s',userVars2.electrodeImp), 'fontsize',fontsize,'fontweight','bold','horizontalalignment','left');

y = y - lineinterval;
uicontrol('parent',hFig,'style','text','units','pixel','position',[x y w 22],'string',sprintf('Recording Depth: %s um',userVars2.recordingDepth), 'fontsize',fontsize,'fontweight','bold','horizontalalignment','left');

y = y - lineinterval;
uicontrol('parent',hFig,'style','text','units','pixel','position',[x y w 22],'string',sprintf('Threshold level: %s mV',userVars2.thresholdSetting), 'fontsize',fontsize,'fontweight','bold','horizontalalignment','left');

y = y - lineinterval*2.5;
uicontrol('parent',hFig,'style','text','units','pixel','position',[x y w-20 22*2.5],'string',sprintf('Activity comments: %s',userVars2.activityComments), 'fontsize',fontsize,'fontweight','bold','horizontalalignment','left');

if isfield(userVars2,'comments')
  y = y - lineinterval*3;
  uicontrol('parent',hFig,'style','text','units','pixel','position',[x y w-20 22*3],'string',sprintf('General Comments: %s',userVars2.comments), 'fontsize',fontsize,'fontweight','bold','horizontalalignment','left');
end

end

load('imageEyeMap_Predata');

if isfield(plotSwitch,'imageEyeMap') && plotSwitch.imageEyeMap
  % Add Colors to each Line
  % Add footage underneath
  % times = -psthPre:psthImDur+psthPost;
    stimStartInd = psthPre+lfpPaddedBy + 1;
    stimEndInd = stimStartInd + psthImDur - 1;
  
  for ii = 1:length(analogInByEvent)
    %Isolate the eye data for a stimuli, its name, and path.
    eyeInByEvent = squeeze(analogInByEvent{ii}(:,1:2,:,stimStartInd:stimEndInd)); %Grab the eye signal for an event
    stimInfo = dir(strcat(stimDir, '/**/', translationTable{ii})); %use the translationTable to find the file
    stimPath = [stimInfo(1).folder filesep stimInfo(1).name]; %create its path.
    screenSize = taskData.screenStats.screenSize;
    PixelsPerDegree = taskData.screenStats.PixelsPerDegree;
    DiagonalSize = taskData.screenStats.DiagonalSize;
    
    %Open the Video, Get some Info on it
    stimVid = VideoReader(stimPath); %#ok<TNMLP> %Open the file.
    TotalFrames = 83;
    samplesPerFrame = floor(psthImDur/TotalFrames);
    stimPos = [-stimVid.Width/2 -stimVid.Height/2];
    
    %Down sample eye signal so that we have one point per frame, convert to
    %pixels
    eyeInByEventDownsampled = zeros(2, size(eyeInByEvent,2), TotalFrames);
    for jj = 1:size(eyeInByEventDownsampled, 2)
      for kk = 1:TotalFrames
        startInd = 1 + ((kk - 1) * samplesPerFrame);
        endInd =  startInd + samplesPerFrame;
        eyeInByEventDownsampled(1,jj,kk) = mean(eyeInByEvent(1,jj,startInd:endInd));
        eyeInByEventDownsampled(2,jj,kk) = mean(eyeInByEvent(2,jj,startInd:endInd));
      end
    end
    eyeInByEventDownsampled(1,:,:) = eyeInByEventDownsampled(1,:,:) * PixelsPerDegree(1);
    eyeInByEventDownsampled(2,:,:) = eyeInByEventDownsampled(2,:,:) * PixelsPerDegree(2);    
    
    
    %Draw each Frame with its Eye signal, from all trials. Cycle through
    %all the Frames, adding the appropriate points to each line. 
    hf = figure;
    title(sprintf('%s eye signal', translationTable{ii}));
    set(hf, 'position', [100 100 screenSize{1} screenSize{2}])
    stimVid.CurrentTime = 0;
    xlim([-screenSize{1}/2 screenSize{1}/2])
    ylim([-screenSize{2}/2 screenSize{2}/2])
    hold on
    axis ij

    %initialize AnimatedLines array
    imagesc(readFrame(stimVid), 'Xdata', stimPos(1), 'YData', stimPos(2));
    A = cell(size(eyeInByEventDownsampled, 2), 1);
    for jj = 1:length(A)
      A{jj} = animatedline('color' ,colors{mod(jj,length(colors))+1}, 'lineWidth', 2 );
    end
    
    %Draw each Frame by drawing from the stimVid, then drawing all the
    %lines.
    for jj = 2:TotalFrames
      imagesc(readFrame(stimVid), 'Xdata', stimPos(1), 'YData', stimPos(2));
      for kk = 1:length(A)
        addpoints(A{kk}, eyeInByEventDownsampled(1,kk,jj), eyeInByEventDownsampled(2,kk,jj))
      end
      hf.Children.Children = circshift(hf.Children.Children,-1);
      pause(1/stimVid.FrameRate)
    end    
  end
end
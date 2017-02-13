%% Converts offline sorter variables into a structure for easy manipulation

for i = 1:96
    ChannelName = ['Chan', num2str(i, '%03.0f')];
    if exist(ChannelName, 'var')
        VariableName = eval(ChannelName);
        Channels.(['Chan', num2str(i, '%3.0f')]) = VariableName;
    end
end

clear Chan0*;
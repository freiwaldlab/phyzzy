function [ ] = buildStimParamFileSocialInt()
% BuildStimParamFileSocialInt
%Creates a .m file containing cell arrays detailing the parameters of a
%particular stimulus set. Julia produced 4 orders of the videos detailed
%below.

%File produced 2 Cell Arrays. The first is a ParamArray, with the name of
%the video. {{VideoID1; {label1}; {label2}...}{VideoID2; {label1};
%{label2}...}. The second is a TriggerLabels, containing all possible labels. 

categoryLabels = {'objects', 'twomonkeys', 'interaction', 'socialInteraction', ...
    'fighting', 'mounting', 'grooming', 'chasing','goalDirected' ...
    'monkeyscramble','scene','scenescramble', 'scramble' };

paramArray = {...
{'Order1Interaction_fighting.avi'; 'twomonkeys'; 'socialInteraction'; 'fighting'}; ...
{'Dephased_Order1Interaction_chasing.avi'; 'monkeyscramble'}; ...
{'Order1Interaction_mounting.avi'; 'twomonkeys'; 'socialInteraction'; 'mounting'}; ...
{'Order1Interaction_grooming2.avi'; 'twomonkeys'; 'socialInteraction'; 'grooming'}; ...
{'Order1MonkLeft_Order1MonkRight_5.avi'; 'twomonkeys'}; ...
{'Order1GoalLeft_Order1GoalRight_2.avi'; 'twomonkeys'; 'goalDirected'}; ...
{'Order2LandscapesFull_3.avi'; 'scene' }; ...
{'Order1MonkLeft_Order1MonkRight_1.avi'; 'twomonkeys'}; ...
{'Dephased_Order1Interaction_fighting.avi'; 'monkeyscramble'}; ...
{'Dephased_Order1Interaction_grooming2.avi'; 'monkeyscramble'}; ...
{'Order1Interaction_chasing.avi'; 'twomonkeys'; 'socialInteraction'; 'chasing'}; ...
{'Order2LandscapesFull_1.avi'; 'scene' }; ...
{'Dephased_Order1Interaction_grooming1.avi'; 'monkeyscramble'}; ...
{'Order1GoalLeft_Order1GoalRight_1.avi'; 'twomonkeys'; 'goalDirected'}; ...
{'Dephased_Order1Interaction_mounting.avi'; 'monkeyscramble'}; ...
{'obj_interact1_2.avi'; 'objects'; 'interaction'}; ...
{'obj_interact1_1.avi'; 'objects'; 'interaction'}; ...
{'Order1Interaction_grooming1.avi'; 'twomonkeys'; 'socialInteraction'; 'grooming'}};

pictureLabels = {...
    'fighting1'; ...
    'scrambleChasing'; ...
    'mounting1'; ...
    'grooming2'; ...
    'nonSocialIdleMonkeys1'; ...
    'nonSocialGoalMonkeys1'; ...
    'landscapesFull3'; ...
    'nonSocialIdleMonkeys2'; ...
    'scrambleFighting1'; ...
    'scrambleGrooming2'; ...
    'chasing1'; ...
    'landscapesFull1'; ...
    'scrambleGrooming1'; ...
    'nonSocialGoalMonkeys2'; ...
    'scrambleMounting1'; ...
    'objInteract2'; ...
    'objInteract1'; ...
    'grooming1'; ...
    };

save('StimParamsFullSocialInt.mat','paramArray','categoryLabels','pictureLabels')



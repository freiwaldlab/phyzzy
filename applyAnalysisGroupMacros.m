function [analysisGroups] = applyAnalysisGroupMacros(analysisGroups, categoryLabels, eventLabels, eventCategories, colors)
%applyAnalysisGroupMacros unpacks the following macros in analysisGroups:
%     - ALL_CATS
%     - ALL_CATS_SPLIT
%     - ALL_EVENTS
%     - ALL_EVENTS_SPLIT
%     - ALL_EVENTS_IN_* (where * is a category label; e.g. ALL_EVENTS_IN_face)
%     - ALL_EVENTS_IN_ALL_CATS_WITH_CAT_AVE
%     - ALL_EVENTS_IN_ALL_CATS
%
%     ALL_CATS, ALL_EVENTS, and ALL_EVENTS_IN_* add items to the group
%       in which they appear. The new items are inserted in place of the
%       macro, so e.g. with categoryLabels = {'face','object','body'},
%       {'face1','ALL_CATS','face2'} becomes {'face1','face','object','body','face2'}
%
%     ALL_CATS_SPLIT and ALL_EVENTS_SPLIT create a new group for each event/category. 
%       Insertion occurs as in ALL_CATS and ALL_EVENTS
%
%     ALL_EVENTS_IN_ALL_CATS_WITH_CAT_AVE and ALL_EVENTS_IN_ALL_CATS make new groups,
%       one group per category. Each group contains the elements of the
%       original group, wtih the new entries inserted in place of the macro
%       as in ALL_CATS etc. 
%     
%     Macros can appear in any combination and number within the group
%       structure. They will be unpacked from left-to-right. So, for
%       instance, {ALL_CATS,ALL_EVENTS_IN_ALL_CATS} would lead to one group per
%       category, and each group would consist of all the categories
%       followed by all the events in one category.
%
%     When a macro appears at the end of a group, it does not need to have
%       an associated color, but if entries (macros or not) follow the marco,
%       if must have an associated dummy entry in the colors cell to maintain indexing.
%       So e.g. groups = {{'face','ALL_EVENTS'}}, colors = {{'r'}} is okay. 
%       And {{'face','ALL_EVENTS','object}}, colors = {{'r','','b'}} is okay. 
%       But {{'face','ALL_EVENTS','object}}, colors = {{'r','b'}} will cause fatal errors.
%
%     Note: macros are not yet implemented for groups with depth > 1. This
%           function will leave those groups unchanged.
%
%     When a macro adds entries to a group, it adds colors via modular
%       indexing into the cell array colors, which can be provided as an
%       argument, and will otherwise default to the Matlab colorspec optons.
%   
%     Parameters:
%       - analysisGroups: struct with arbitrary fields, each of which is a struct
%                         with (at least) fields groups, colors, and names.
%                         Groups and colors are cell arrays of cell arrays.
%                         Names is a cell array. 
%       - categoryLabels: cell array of strings, optional if only macro
%                         present is ALL_EVENTS
%       - eventLabels: cell array of strings, optional if onlt macro
%                      present is ALL_CATS
%       - eventCategories: cell array of cell arrays as strings, indexed to
%                          match eventLabels. Optional if only macros
%                          present are ALL_EVENTS and ALL_CATS
%       - colors: optional, cell array of either colorspec letters or numeric rgb triplets
%
%
%     Possible extensions:
%       Add a parameter that specifies colors for certain, or all,
%       categories and events
%       
%       Add handling for depth = 2 groups


if ~exist('colors','var')
  colors = {'b';'c';'y';'g';'m';'r';'k'};
end

analysisGroupsFields = fieldnames(analysisGroups);
for field_i = 1:length(analysisGroupsFields)
  fieldName = analysisGroupsFields{field_i};
  field = analysisGroups.(fieldName);
  group_i = 0;
  while group_i < length(field.groups)
    group_i = group_i + 1;
    item_i = 0;
    while item_i < length(field.groups{group_i})
      item_i = item_i + 1;
      item = field.groups{group_i}{item_i};
      if iscell(item)  % skip depth > 1 groups; handling not implemented
        continue
      end
      if strcmp(item,'ALL_CATS')
        newGroup = vertcat(field.groups{group_i}(1:item_i-1),reshape(categoryLabels,length(categoryLabels),1));
        newColors = vertcat(field.colors{group_i}(1:item_i-1),colors(mod(0:length(categoryLabels)-1,length(colors)-1)+1));
        if item_i < length(field.groups{group_i})
          newGroup = vertcat(newGroup,field.groups{group_i}(item_i+1:end));
          newColors = vertcat(newColors,field.colors{group_i}(item_i+1:end));
        end
        field.groups{group_i} = newGroup;
        field.colors{group_i} = newColors;

      elseif strcmp(item,'ALL_EVENTS')
        newGroup = vertcat(field.groups{group_i}(1:item_i-1),reshape(eventLabels,length(eventLabels),1));
        newColors = vertcat(field.colors{group_i}(1:item_i-1),colors(mod(0:length(eventLabels)-1,length(colors)-1)+1));
        if item_i < length(field.groups{group_i})
          newGroup = vertcat(newGroup,field.groups{group_i}(item_i+1:end));
          newColors = vertcat(newColors,field.colors{group_i}(item_i+1:end));
        end
        field.groups{group_i} = newGroup;
        field.colors{group_i} = newColors;

      elseif ~isempty(regexp(item,regexptranslate('wildcard','ALL_EVENTS_IN_*'), 'ONCE'))
        catLabel = item(15:end);
        newGroup = field.groups{group_i}(1:item_i-1);
        newColors = field.colors{group_i}(1:item_i-1);
        eventsAdded = 0;
        for event_i = 1:length(eventLabels)
          if any(strcmp(eventCategories{event_i},catLabel))
            newGroup = vertcat(newGroup,eventLabels{event_i});
            eventsAdded = eventsAdded + 1;
          end
        end
        newColors = vertcat(newColors,colors(mod(0:eventsAdded-1,length(colors)-1)+1));
        if item_i < length(field.groups{group_i})
          newGroup = vertcat(newGroup,field.groups{group_i}(item_i+1:end));
          newColors = vertcat(newColors,field.colors{group_i}(item_i+1:end));
        end
        field.groups{group_i} = newGroup;
        field.colors{group_i} = newColors;

      elseif strcmp(item,'ALL_EVENTS_IN_ALL_CATS_WITH_CAT_AVE')         
        newGroupPrefix = field.groups{group_i}(1:item_i-1);
        newColorPrefix = field.colors{group_i}(1:item_i-1);
        if item_i < length(field.groups{group_i})
          newGroupSuffix = field.groups{group_i}(item_i+1:end);
          newColorSuffix = field.colors{group_i}(item_i+1:end);
        else
          newGroupSuffix = {};
          newColorSuffix = {};
        end
        newGroups = field.groups(1:group_i-1);
        newColors = field.colors(1:group_i-1);
        newNames = field.names(1:group_i-1);
        if group_i < length(field.groups)
          newGroupsSuffix = vertcat(newGroups,field.groups(group_i+1:end));
          newColorsSuffix = vertcat(newColors,field.colors(group_i+1:end));
          newNamesSuffix = vertcat(newNames,field.names(group_i+1:end));
        else
          newGroupsSuffix = {};
          newColorsSuffix = {};
          newNamesSuffix = {};
        end
        for cat_i = 1:length(categoryLabels)
          catLabel = categoryLabels{cat_i};
          newGroup = vertcat(newGroupPrefix,catLabel);
          eventsAdded = 1;
          for event_i = 1:length(eventLabels)
            if any(strcmp(eventCategories{event_i},catLabel))
              eventsAdded = eventsAdded + 1;
              newGroup = vertcat(newGroup,eventLabels{event_i});
            end
          end
          newGroup = vertcat(newGroup, newGroupSuffix);
          newGroups = vertcat(newGroups, {newGroup});
          newColor = vertcat(newColorPrefix,colors(mod(0:eventsAdded-1,length(colors)-1)+1),newColorSuffix);
          newColors = vertcat(newColors,{newColor});
          newNames = vertcat(newNames, strcat(catLabel,'s'));
        end
        newGroups = vertcat(newGroups, newGroupsSuffix);
        newColors = vertcat(newColors, newColorsSuffix);
        newNames = vertcat(newNames, newNamesSuffix);
        field.groups = newGroups;
        field.colors = newColors;
        field.names = newNames;
     
      elseif strcmp(item,'ALL_EVENTS_IN_ALL_CATS')
        newGroupPrefix = field.groups{group_i}(1:item_i-1);
        newColorPrefix = field.colors{group_i}(1:item_i-1);
        if item_i < length(field.groups{group_i})
          newGroupSuffix = field.groups{group_i}(item_i+1:end);
          newColorSuffix = field.colors{group_i}(item_i+1:end);
        else
          newGroupSuffix = {};
          newColorSuffix = {};
        end
        newGroups = field.groups(1:group_i-1);
        newColors = field.colors(1:group_i-1);
        newNames = field.names(1:group_i-1);
        if group_i < length(field.groups)
          newGroupsSuffix = vertcat(newGroups,field.groups(group_i+1:end));
          newColorsSuffix = vertcat(newColors,field.colors(group_i+1:end));
          newNamesSuffix = vertcat(newNames,field.names(group_i+1:end));
        else
          newGroupsSuffix = {};
          newColorsSuffix = {};
          newNamesSuffix = {};
        end
        for cat_i = 1:length(categoryLabels)
          catLabel = categoryLabels{cat_i};
          newGroup = newGroupPrefix;
          eventsAdded = 0;
          for event_i = 1:length(eventLabels)
            if any(strcmp(eventCategories{event_i},catLabel))
              eventsAdded = eventsAdded + 1;
              newGroup = vertcat(newGroup,eventLabels{event_i});
           end
          end
          newGroup = vertcat(newGroup, newGroupSuffix);
          newGroups = vertcat(newGroups, {newGroup});
          newColor = vertcat(newColorPrefix,colors(mod(0:eventsAdded-1,length(colors)-1)+1),newColorSuffix);
          newColors = vertcat(newColors,{newColor});
          newNames = vertcat(newNames, strcat(catLabel,'s'));
        end
        newGroups = vertcat(newGroups, newGroupsSuffix);
        newColors = vertcat(newColors, newColorsSuffix);
        field.groups = newGroups; 
        field.colors = newColors;
        field.names = newNames;
        
      elseif strcmp(item,'ALL_CATS_SPLIT')
        newGroupPrefix = field.groups{group_i}(1:item_i-1);
        newColorPrefix = field.colors{group_i}(1:item_i-1);
        if item_i < length(field.groups{group_i})
          newGroupSuffix = field.groups{group_i}(item_i+1:end);
          newColorSuffix = field.colors{group_i}(item_i+1:end);
        else
          newGroupSuffix = {};
          newColorSuffix = {};
        end
        newGroups = field.groups(1:group_i-1);
        newColors = field.colors(1:group_i-1);
        newNames = field.names(1:group_i-1);
        if group_i < length(field.groups)
          newGroupsSuffix = vertcat(newGroups,field.groups(group_i+1:end));
          newColorsSuffix = vertcat(newColors,field.colors(group_i+1:end));
          newNamesSuffix = vertcat(newNames,field.names(group_i+1:end));
        else
          newGroupsSuffix = {};
          newColorsSuffix = {};
          newNamesSuffix = {};
        end
        for cat_i = 1:length(categoryLabels)
          catLabel = categoryLabels{cat_i};
          newGroup = vertcat(newGroupPrefix,catLabel);
          newGroup = vertcat(newGroup, newGroupSuffix);
          newGroups = vertcat(newGroups, {newGroup});
          newColor = vertcat(newColorPrefix,colors(1),newColorSuffix);
          newColors = vertcat(newColors,{newColor});
          newNames = vertcat(newNames, strcat(catLabel));
        end
        newGroups = vertcat(newGroups, newGroupsSuffix);
        newColors = vertcat(newColors, newColorsSuffix);
        newNames = vertcat(newNames, newNamesSuffix);
        field.groups = newGroups;
        field.colors = newColors;
        field.names = newNames;
        
      elseif strcmp(item,'ALL_EVENTS_SPLIT')
        newGroupPrefix = field.groups{group_i}(1:item_i-1);
        newColorPrefix = field.colors{group_i}(1:item_i-1);
        if item_i < length(field.groups{group_i})
          newGroupSuffix = field.groups{group_i}(item_i+1:end);
          newColorSuffix = field.colors{group_i}(item_i+1:end);
        else
          newGroupSuffix = {};
          newColorSuffix = {};
        end
        newGroups = field.groups(1:group_i-1);
        newColors = field.colors(1:group_i-1);
        newNames = field.names(1:group_i-1);
        if group_i < length(field.groups)
          newGroupsSuffix = vertcat(newGroups,field.groups(group_i+1:end));
          newColorsSuffix = vertcat(newColors,field.colors(group_i+1:end));
          newNamesSuffix = vertcat(newNames,field.names(group_i+1:end));
        else
          newGroupsSuffix = {};
          newColorsSuffix = {};
          newNamesSuffix = {};
        end
        for event_i = 1:length(eventLabels)
          eventLabel = eventLabels{event_i};
          newGroup = vertcat(newGroupPrefix,eventLabel);
          newGroup = vertcat(newGroup, newGroupSuffix);
          newGroups = vertcat(newGroups, {newGroup});
          newColor = vertcat(newColorPrefix,colors(1),newColorSuffix);
          newColors = vertcat(newColors,{newColor});
          newNames = vertcat(newNames, strcat(eventLabel));
        end
        newGroups = vertcat(newGroups, newGroupsSuffix);
        newColors = vertcat(newColors, newColorsSuffix);
        newNames = vertcat(newNames, newNamesSuffix);
        field.groups = newGroups;
        field.colors = newColors;
        field.names = newNames;
      else
        continue
      end          
    end
  end
  analysisGroups.(fieldName) = field;
end
end


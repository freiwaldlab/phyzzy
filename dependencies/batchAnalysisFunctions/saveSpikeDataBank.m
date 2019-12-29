function output = saveSpikeDataBank(structData, parts, SaveOrLoad, outDir)
% Save or loads large structs into slices. used by preprocessBatchAnalysis
% and batchAnalysis for intermittent saving.
% - structData: The large struct
% - parts: the number of parts to divide it into.
% - SaveOrLoad: whether to 'save' or 'load' the struct.
switch SaveOrLoad
  case 'save'
    %saves variable in N parts.
    fieldsToStore = fields(structData);
    fieldCount = length(fields(structData));
    fieldsPerFile = ceil(fieldCount/parts);
    output = cell(ceil(fieldCount/fieldsPerFile),1);
    spikeSliceCount = 1;
    spikeDataBankSlice = struct();
    
    for field_i = 1:length(fieldsToStore)
      spikeDataBankSlice.(fieldsToStore{field_i}) = structData.(fieldsToStore{field_i});
      if mod(field_i, fieldsPerFile) == 0 || field_i == length(fieldsToStore)
        output{spikeSliceCount} = fullfile(outDir, sprintf('spikeDataBank_%d.mat', spikeSliceCount));
        fprintf('Saving %s now...', output{spikeSliceCount});
        save(output{spikeSliceCount},'spikeDataBankSlice')
        clear spikeDataBankSlice
        spikeSliceCount = spikeSliceCount + 1;
      end
    end
  case 'load'
    spikeDataBankTmp = struct();
    filesToCombine = dir(fullfile(outDir, ['spikeDataBank_*.mat']));
    for ii = 1:length(filesToCombine)
      tmp = load(filesToCombine(ii).name);
      fieldsToCombine = fields(tmp.spikeDataBankSlice);
      for field_ind = 1:length(fieldsToCombine)
        spikeDataBankTmp.(fieldsToCombine{field_ind}) = tmp.spikeDataBankSlice.(fieldsToCombine{field_ind});
      end
    end
    output = spikeDataBankTmp;
end
end
%% prfSampleModel

% addpath(genpath('/home/hanseohe/Documents/GitHub/stimulusVignetting'))
% for sub = 7:8
%     fprintf('%s. %d. %d ...\n','prfSampleModel',sub,1);
%     prfSampleModel(sub,1);
%     fprintf('%s. %d. %d ...\n','prfSampleModel',sub,2);
%     prfSampleModel(sub,2);
%     fprintf('%s. %d. %d ...\n','prfSampleModel',sub,3);
%     prfSampleModel(sub,3);
%     fprintf('%s. %d. %d ...\n','prfSampleModel',sub,4);
%     prfSampleModel(sub,4);
% end

%% regressPrfSplit
for sub = 2:8
fprintf('%s. %d. %d ...\n','regressPrfSplit',sub, 1);
regressPrfSplit_sfmean(sub, 1);
fprintf('%s. %d. %d ...\n','regressPrfSplit',sub, 2);
regressPrfSplit_sfmean(sub, 2);
fprintf('%s. %d. %d ...\n','regressPrfSplit',sub, 3);
regressPrfSplit_sfmean(sub, 3);
fprintf('%s. %d. %d ...\n','regressPrfSplit',sub, 4);
regressPrfSplit_sfmean(sub, 4);

end

%% getVoxPref
% for sub = 1:8
%     getVoxPref_regress(sub,4, 3);
%     fprintf('%s. %d. %d ...\n','getVoxPref_regress',sub, 3);
% end
%% prfSampleModel

addpath(genpath('/home/hanseohe/Documents/GitHub/stimulusVignetting'))
for sub = 3:8
    fprintf('%s. %d. %d ...\n','prfSampleModel',sub,1);
    prfSampleModel(sub,1);
    fprintf('%s. %d. %d ...\n','prfSampleModel',sub,2);
    prfSampleModel(sub,2);
    fprintf('%s. %d. %d ...\n','prfSampleModel',sub,3);
    prfSampleModel(sub,3);
    fprintf('%s. %d. %d ...\n','prfSampleModel',sub,4);
    prfSampleModel(sub,4);
end

% %% regressPrfSplit
% for sub = 5:8
% fprintf('%s. %d. %d ...\n','regressPrfSplit',sub, 1);
% regressPrfSplit(sub, 1);
% fprintf('%s. %d. %d ...\n','regressPrfSplit',sub, 2);
% regressPrfSplit(sub, 2);
% fprintf('%s. %d. %d ...\n','regressPrfSplit',sub, 3);
% regressPrfSplit(sub, 3);
% fprintf('%s. %d. %d ...\n','regressPrfSplit',sub, 4);
% regressPrfSplit(sub, 4);
% 
% end
% 
% %% 
% getVoxPref(1,4)
% getVoxPref(2,4)
% getVoxPref(3,4)
% getVoxPref(4,4)
% getVoxPref(5,4); 
% getVoxPref(6,4); 
% getVoxPref(7,4); 
% getVoxPref(8,4);
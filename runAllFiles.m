%% prfSampleModel
for sub = 2:8
    fprintf('%s. %d. %d ...\n','prfSampleModel',sub,1);
    prfSampleModel(sub,1);
    fprintf('%s. %d. %d ...\n','prfSampleModel',sub,2);
    prfSampleModel(sub,2);
    fprintf('%s. %d. %d ...\n','prfSampleModel',sub,3);
    prfSampleModel(sub,3);
    fprintf('%s. %d. %d ...\n','prfSampleModel',sub,4);
    prfSampleModel(sub,4);
end

%% regressPrfSplit
% for sub = 2:8
% 
% fprintf('%s. %d. %d ...\n','regressPrfSplit',sub,1);
% regressPrfSplit(sub, 4);
% 
% end
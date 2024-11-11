%% prfSampleModel
% for sub = 1:1
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
for sub = 1:8
fprintf('%s. %d. %d ...\n','regressPrfSplit_control',sub, 1);
regressPrfSplit_control(sub, 1);
fprintf('%s. %d. %d ...\n','regressPrfSplit_control',sub, 2);
regressPrfSplit_control(sub, 2);
fprintf('%s. %d. %d ...\n','regressPrfSplit_control',sub, 3);
regressPrfSplit_control(sub, 3);
fprintf('%s. %d. %d ...\n','regressPrfSplit_control',sub, 4);
regressPrfSplit_control(sub, 4);
end

%% getVoxPref_control
% for sub = 2:8
%     getVoxPref_control(sub,4);
% end
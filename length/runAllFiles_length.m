%% prfSampleModel_length

addpath(genpath('/home/hanseohe/Documents/GitHub/stimulusVignetting'))
for sub = 1:1
    tic;
    fprintf('%s. %d. %d ...\n','prfSampleModel',sub,1);
    prfSampleModel_length(sub,1);
    toc;

    tic;
    fprintf('%s. %d. %d ...\n','prfSampleModel',sub,2);
    prfSampleModel_length(sub,2);
    toc;

    tic;
    fprintf('%s. %d. %d ...\n','prfSampleModel',sub,3);
    prfSampleModel_length(sub,3);
    toc;

    tic;
    fprintf('%s. %d. %d ...\n','prfSampleModel',sub,4);
    prfSampleModel_length(sub,4);
    toc;
end

%% regressPrfSplit
% for sub = 5:8
% fprintf('%s. %d. %d ...\n','regressPrfSplit',sub, 1);
% regressPrfSplit_new(sub, 1);
% fprintf('%s. %d. %d ...\n','regressPrfSplit',sub, 2);
% regressPrfSplit_new(sub, 2);
% fprintf('%s. %d. %d ...\n','regressPrfSplit',sub, 3);
% regressPrfSplit_new(sub, 3);
% fprintf('%s. %d. %d ...\n','regressPrfSplit',sub, 4);
% regressPrfSplit_new(sub, 4);
% 
% end
% 
% %% 
% getVoxPref_new(1,4)
% getVoxPref_new(2,4)
% getVoxPref_new(3,4)
% getVoxPref_new(4,4)
% getVoxPref_new(5,4); getVoxPref_new(6,4); getVoxPref_new(7,4); getVoxPref_new(8,4);
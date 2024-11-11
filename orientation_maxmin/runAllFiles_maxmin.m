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
for type = 1:2

        % for sub = 1:8
        %     fprintf('%s. %d. %d. %d ...\n','regressPrfSplit',sub, 1, type);
        %     regressPrfSplit_maxmin(sub, 1, type);
        %     fprintf('%s. %d. %d. %d ...\n','regressPrfSplit',sub, 2, type);
        %     regressPrfSplit_maxmin(sub, 2, type);
        %     fprintf('%s. %d. %d. %d ...\n','regressPrfSplit',sub, 3, type);
        %     regressPrfSplit_maxmin(sub, 3, type);
        %     fprintf('%s. %d. %d. %d ...\n','regressPrfSplit',sub, 4, type);
        %     regressPrfSplit_maxmin(sub, 4, type);
        % 
        % end

    for sub = 1:8
        fprintf('%s. %d. %d. %d ...\n','regressPrfSplit_new',sub, 1, type);
        regressPrfSplit_maxmin_new(sub, 1, type);
        fprintf('%s. %d. %d. %d ...\n','regressPrfSplit_new',sub, 2, type);
        regressPrfSplit_maxmin_new(sub, 2, type);
        fprintf('%s. %d. %d. %d ...\n','regressPrfSplit_new',sub, 3, type);
        regressPrfSplit_maxmin_new(sub, 3, type);
        fprintf('%s. %d. %d. %d ...\n','regressPrfSplit_new',sub, 4, type);
        regressPrfSplit_maxmin_new(sub, 4, type);

    end

    for sub = 1:8
        fprintf('%s. %d. %d. %d ...\n','regressPrfSplit_control',sub, 1, type);
        regressPrfSplit_maxmin_control(sub, 1, type);
        fprintf('%s. %d. %d. %d ...\n','regressPrfSplit_control',sub, 2, type);
        regressPrfSplit_maxmin_control(sub, 2, type);
        fprintf('%s. %d. %d. %d ...\n','regressPrfSplit_control',sub, 3, type);
        regressPrfSplit_maxmin_control(sub, 3, type);
        fprintf('%s. %d. %d. %d ...\n','regressPrfSplit_control',sub, 4, type);
        regressPrfSplit_maxmin_control(sub, 4, type);
    end
end
%
%% getVoxPref
% for sub = 1:8
%     getVoxPref(sub,4)
% end
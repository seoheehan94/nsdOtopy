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
isub_values = 1:8;
visualRegions_values = 1:4;
pairTypes = {'old_control', 'old_ori', 'control_ori'};
imgTypes = {'Top', 'Bottom'};

% Initialize regressPrfSplit_maxminPatch as needed, e.g., as a struct or cell array

% Loop through all combinations of isub, visualRegions, pairType, and imgType
for isub = isub_values
    for visualRegion = visualRegions_values
        for pairTypeIdx = 1:length(pairTypes)
            for imgTypeIdx = 1:length(imgTypes)
                
                pairType = pairTypes{pairTypeIdx};
                imgType = imgTypes{imgTypeIdx};

                fprintf('%s. %d. %d. %s. %s. ...\n','regressPrfSplit_maxminPatch',isub,visualRegion, pairType, imgType);
                regressPrfSplit_maxminPatch(isub,visualRegion,pairType, imgType)
                
            end
        end
    end
end
%
%% getVoxPref
% for sub = 1:8
%     getVoxPref(sub,4)
% end
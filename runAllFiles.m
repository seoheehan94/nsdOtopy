%% prfSampleModel_new
for sub = 1:4
    fprintf('%s. %d. %d ...\n','prfSampleModel',sub,1);
    prfSampleModel_new(sub,1);
    fprintf('%s. %d. %d ...\n','prfSampleModel',sub,2);
    prfSampleModel_new(sub,2);
    fprintf('%s. %d. %d ...\n','prfSampleModel',sub,3);
    prfSampleModel_new(sub,3);
    fprintf('%s. %d. %d ...\n','prfSampleModel',sub,4);
    prfSampleModel_new(sub,4);
end

%% regressPrfSplit
for sub = 1:4
fprintf('%s. %d. %d ...\n','regressPrfSplit',sub, 1);
regressPrfSplit_new(sub, 1);
fprintf('%s. %d. %d ...\n','regressPrfSplit',sub, 2);
regressPrfSplit_new(sub, 2);
fprintf('%s. %d. %d ...\n','regressPrfSplit',sub, 3);
regressPrfSplit_new(sub, 3);
fprintf('%s. %d. %d ...\n','regressPrfSplit',sub, 4);
regressPrfSplit_new(sub, 4);

end
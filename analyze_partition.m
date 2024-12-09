clear all;


savefolder = '/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Orientation/analyses/partition/';


%% 1.  mean R2, AIC, BIC %%
totalR2_old = cell(4,1);
totalR2_unique_ori = cell(4,1);
totalR2_combined = cell(4,1);
totalsharedVariance = cell(4,1);

for visualRegion =1:4
    for isub = 1:8
        fprintf('isub:%d. V:%d. ...\n',isub,visualRegion);
        load([savefolder 'regressPartition_oldori_v' num2str(visualRegion) '_sub' num2str(isub) '.mat']);
        if length(rois)>1
            oldNsd = nsd;
            nsd.R2_model1 = [];
            nsd.R2_model2_orth = [];
            nsd.R2_combined = [];
            nsd.sharedVariance = [];
            for iroi=1:length(rois)
                nsd.R2_model1 = cat(2,nsd.R2_model1,oldNsd.R2_model1{iroi});
                nsd.R2_model2_orth = cat(2,nsd.R2_model2_orth,oldNsd.R2_model2_orth{iroi});
                nsd.R2_combined = cat(2,nsd.R2_combined,oldNsd.R2_combined{iroi});
                nsd.sharedVariance = cat(2,nsd.sharedVariance,oldNsd.sharedVariance{iroi});
            end
        end
        nsd.R2_model1(3,:,:) = mean(nsd.R2_model1,1);
        nsd.R2_model2_orth(3,:,:) = mean(nsd.R2_model2_orth,1);
        nsd.R2_combined(3,:,:) = mean(nsd.R2_combined,1);
        nsd.sharedVariance(3,:,:) = mean(nsd.sharedVariance,1);
        


        %% total values
        totalR2_old{visualRegion}{end+1} = nsd.R2_model1(3,:,:);
        totalR2_unique_ori{visualRegion}{end+1} = nsd.R2_model2_orth(3,:,:);
        totalR2_combined{visualRegion}{end+1} = nsd.R2_combined(3,:,:);
        totalsharedVariance{visualRegion}{end+1} = nsd.sharedVariance(3,:,:);
    end

end


%% average participants
allroiR2_old=[];
V1R2_old=[];
allroiR2_unique_ori=[];
V1R2_unique_ori=[];
allroiR2_combined=[];
V1R2_combined=[];
allroisharedVariance=[];
V1sharedVariance=[];
for visualRegion =1:4
    curRoiR2_old = [];
    curV1R2_old = [];
    curRoiR2_unique_ori = [];
    curV1R2_unique_ori = [];
    curRoiR2_combined = [];
    curV1R2_combined = [];
    curRoisharedVariance = [];
    curV1sharedVariance = [];
    
    for isub = 1:8
            curRoiR2_old = [curRoiR2_old, totalR2_old{visualRegion}{isub}];

        curV1R2_old = [curV1R2_old, totalR2_old.(fieldsCon{i}){j}{1}(3,:)];
    end
    % writematrix(curRoiR2OriSplit', ['/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Orientation/analyses/MaxMin/allroiR2', imgType{curimgtype}, '_', (fieldsCon{i}), '.csv']);
    % writematrix(curV1R2OriSplit', ['/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Orientation/analyses/MaxMin/V1R2', imgType{curimgtype}, '_', (fieldsCon{i}), '.csv']);

    allroiR2_old(i) = mean(curRoiR2_old, 'omitnan');
    V1R2_old(i) = mean(curV1R2_old,'omitnan');
end
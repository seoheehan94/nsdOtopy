clear all;


savefolder = '/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Orientation/analyses/partition/';


%% 1.  mean R2, AIC, BIC %%
totalsharedVariance = cell(4,1);
totaluniqueModel_old = cell(4,1);
totaluniqueModel_ori = cell(4,1);

for visualRegion =1:4
    for isub = 1:8
        fprintf('isub:%d. V:%d. ...\n',isub,visualRegion);
        load([savefolder 'regressPartition_oldori_v' num2str(visualRegion) '_sub' num2str(isub) '.mat']);
        if length(rois)>1
            oldNsd = nsd;
            nsd.sharedVariance = [];
            nsd.uniqueModel1 = [];
            nsd.uniqueModel2 = [];
            for iroi=1:length(rois)
                nsd.sharedVariance = cat(2,nsd.sharedVariance,oldNsd.sharedVariance{iroi});
                nsd.uniqueModel1 = cat(2,nsd.uniqueModel1,oldNsd.uniqueModel1{iroi});
                nsd.uniqueModel2 = cat(2,nsd.uniqueModel2,oldNsd.uniqueModel2{iroi});
            end
        end
        nsd.sharedVariance(3,:,:) = mean(nsd.sharedVariance,1);
        nsd.uniqueModel1(3,:,:) = mean(nsd.uniqueModel1,1);
        nsd.uniqueModel2(3,:,:) = mean(nsd.uniqueModel2,1);
        


        %% total values
        totalsharedVariance{visualRegion}{end+1} = nsd.sharedVariance(3,:,:);
        totaluniqueModel_old{visualRegion}{end+1} = nsd.uniqueModel1(3,:,:);
        totaluniqueModel_ori{visualRegion}{end+1} = nsd.uniqueModel2(3,:,:);
    end

end


%% average participants
fieldsCon = fieldnames(totalsharedVariance);
allroiR2OriSplit=[];
V1R2OriSplit=[];
for i = 1:numel(fieldsCon)
    curRoiR2OriSplit = [];
    curV1R2OriSplit = [];
    
    for j = 1:numel(totalsharedVariance.(fieldsCon{i}))
        for k = 1:size(roiNsdOriR2,2)
            curRoiR2OriSplit = [curRoiR2OriSplit, totalsharedVariance.(fieldsCon{i}){j}{k}(3,:)];
        end

        curV1R2OriSplit = [curV1R2OriSplit, totalsharedVariance.(fieldsCon{i}){j}{1}(3,:)];
    end
    % writematrix(curRoiR2OriSplit', ['/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Orientation/analyses/MaxMin/allroiR2', imgType{curimgtype}, '_', (fieldsCon{i}), '.csv']);
    % writematrix(curV1R2OriSplit', ['/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Orientation/analyses/MaxMin/V1R2', imgType{curimgtype}, '_', (fieldsCon{i}), '.csv']);

    allroiR2OriSplit(i) = mean(curRoiR2OriSplit, 'omitnan');
    V1R2OriSplit(i) = mean(curV1R2OriSplit,'omitnan');
end
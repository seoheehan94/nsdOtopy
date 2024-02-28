%getAverage_length.m

%   uses files created by: regressPrfSplit_length.m
%   creates files used by:
for isub = 2:8
    clearvars -except isub
    %% set up
    prffolder = ['/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Length/prfsample_Len/'];

    bandpass = 1; bandMin = 1; bandMax = 1;
    bandpassStr = '';
    if bandpass
        bandpassStr = ['_bandpass' num2str(bandMin) 'to' num2str(bandMax)];
    end

    mean_split1 = struct;
    mean_split2 = struct;
    mean_all = struct;
    allLenCoef = struct;
    maxCoefLen = struct;

    %% load file
    for visualRegion = 1:4

        load(fullfile(prffolder,['regressPrfSplit' bandpassStr '_v' num2str(visualRegion) '_sub' num2str(isub)  '.mat']), ...
            'nsd', 'numLevels', 'numLengths','rois','nvox','roiPrf','nsplits');

        if visualRegion == 4

            mean_split1.v4 = squeeze(mean(nsd.voxLenCoef{1}(1,:,1:8),2));
            mean_split2.v4 = squeeze(mean(nsd.voxLenCoef{1}(2,:,1:8),2));

            % average splits    
            nsd.voxLenCoef{1}(nsplits+1,:,:) = mean(nsd.voxLenCoef{1},1);
            mean_all.v4 = squeeze(mean(nsd.voxLenCoef{1}(3,:,1:8),2));
            allLenCoef.v4 = nsd.voxLenCoef{1}(:,:,1:8);

        else %for other regions, combine ventral and dorsal
            % combine ventral and dorsal    
            oldNsd = nsd;
            nsd.voxLenCoef{1} = [];
            nsd.voxLenCoef{2} = [];
            for iroi=1:length(rois)
                nsd.voxLenCoef{1} = cat(2,nsd.voxLenCoef{1},oldNsd.voxLenCoef{iroi});
            end

            thisfield = ['v', num2str(visualRegion)];
            mean_split1.(thisfield) = squeeze(mean(nsd.voxLenCoef{1}(1,:,1:8),2));
            mean_split2.(thisfield) = squeeze(mean(nsd.voxLenCoef{1}(2,:,1:8),2));

            % average splits    
            nsd.voxLenCoef{1}(nsplits+1,:,:) = mean(nsd.voxLenCoef{1},1);
            mean_all.(thisfield) = squeeze(mean(nsd.voxLenCoef{1}(3,:,1:8),2));
            allLenCoef.(thisfield) = nsd.voxLenCoef{1}(:,:,1:8);
        end


    end
    %% check consistency between splits
    % for visualRegion = 1:4
    %     thisfield = ['v', num2str(visualRegion)];
    % 
    %     figure;
    %     subplot(1,2,1)
    %     bar(mean_split1.(thisfield),1);
    %     title('split1');
    % 
    %     subplot(1,2,2)
    %     bar(mean_split2.(thisfield),1);
    %     title('split2');
    % 
    %     totalTitle = ['sub', num2str(isub), ' V', num2str(visualRegion)];
    %     sgtitle(totalTitle);
    % 
    %     saveas(gcf,['/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Length/len_hist/split_sub', num2str(isub), 'v', num2str(visualRegion), '.png']);
    % end

    %% get average for each ROI
    % figure;
    % subplot(2,2,1)
    % bar(mean_all.v1,1);
    % title('v1');
    % 
    % subplot(2,2,2)
    % bar(mean_all.v2,1);
    % title('v2');
    % 
    % subplot(2,2,3)
    % bar(mean_all.v3,1);
    % title('v3');
    % 
    % subplot(2,2,4)
    % bar(mean_all.v4,1);
    % title('v4');
    % 
    % totalTitle = ['sub', num2str(isub), ' average for ROI'];
    % sgtitle(totalTitle);

    % saveas(gcf,['/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Length/len_hist/all_sub', num2str(isub), '.png']);


%% get preferred length for each voxel
% save all ROIs to create overlay
roifolder = ['/bwdata/NSDData/nsddata/ppdata/subj0' num2str(isub) '/func1pt8mm/'];
visualRoisFile = fullfile(roifolder,'roi/prf-visualrois.nii.gz');%V1v, V1d, V2v, V2d, V3v, V3d, and hV4
visRoiData = niftiread(visualRoisFile);
roiNames = {'V1v','V1d','V2v','V2d','V3v','V3d','hV4'};
combinedRoiNames = {'V1','V2','V3','hV4'};

ourBrain = visRoiData;
ourBrain(ourBrain == 2) = 1;
ourBrain(ourBrain == 3 | ourBrain == 4) = 2;
ourBrain(ourBrain == 5 | ourBrain == 6) = 3;
ourBrain(ourBrain == 7) = 4;

% find max coef for each voxel
for visualRegion = 1:4
    thisfield = ['v', num2str(visualRegion)];
    [~, nvox, ~] = size(allLenCoef.(thisfield));
    for ivox = 1: nvox
        [~, ilength] = max(allLenCoef.(thisfield)(3,ivox,:));
        maxCoefLen.(thisfield)(1,ivox) = ilength;
    end
end

newBrain = ourBrain;
newBrain(newBrain > 0) = 0;
newBrain(ourBrain == 1) = maxCoefLen.v1(1,:);
newBrain(ourBrain == 2) = maxCoefLen.v2(1,:);
newBrain(ourBrain == 3) = maxCoefLen.v3(1,:);
newBrain(ourBrain == 4) = maxCoefLen.v4(1,:);
newBrain(newBrain < 0) = -1;

for k = 1:4
    curNewBrain = ourBrain;
    curNewBrain(ourBrain ~= k) = -1;
    thisfield = ['v', num2str(k)];
    curNewBrain(ourBrain == k) = maxCoefLen.(thisfield)(1,:);

    newBrainbyROI(:,:,:,k) =curNewBrain;
end

% saveName = ['/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Length/brainVolume/lengthBrain_sub', num2str(isub), '.mat'];
% save(saveName, 'newBrain');

saveName = ['/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Length/brainVolume/lengthBrainbyROI_sub', num2str(isub), '.mat'];
save(saveName, 'newBrainbyROI');

%% save afni file
% size(newBrain)
cd('/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Length/brainVolume');
% %3dinfo betas_session01.nii.gz 
% %3dcalc -a betas_session01.nii.gz[1] -expr a -prefix oneBeta
% command = ['3dinfo betas_session01_sub', num2str(isub), '.nii.gz'];
% system(command);
% command = ['3dcalc -a betas_session01_sub', num2str(isub), '.nii.gz[1] -expr a -prefix oneBeta_sub', num2str(isub)];
% system(command);
% 
currOneBeta = ['oneBeta_sub', num2str(isub), '+orig'];
[err,V,Info] = BrikLoad(currOneBeta);
Info.RootName = ['lengthBrainbyROI_sub', num2str(isub), '+orig'];
opt.Prefix = ['lengthBrainbyROI_sub', num2str(isub)];
WriteBrik(newBrainbyROI,Info,opt);

%% plot figure
% figure;
% subplot(2,2,1)
% histogram(maxCoefLen.v1);
% title('v1');
% 
% subplot(2,2,2)
% histogram(maxCoefLen.v2);
% title('v2');
% 
% subplot(2,2,3)
% histogram(maxCoefLen.v3);
% title('v3');
% 
% subplot(2,2,4)
% histogram(maxCoefLen.v4);
% title('v4');
% 
% totalTitle = ['sub', num2str(isub), ' number of preferred length bin'];
% sgtitle(totalTitle);
% 
% saveas(gcf,['/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Length/len_hist/prefLength_sub', num2str(isub), '.png']);

end
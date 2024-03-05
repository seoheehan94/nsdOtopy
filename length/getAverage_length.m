%getAverage_length.m

%   uses files created by: regressPrfSplit_length.m
%   creates files used by:

roiNames = {'V1v','V1d','V2v','V2d','V3v','V3d','hV4','OPA','PPA','RSC'};
combinedRoiNames = {'V1','V2','V3','hV4','OPA','PPA','RSC'};

for isub = 2:8
     fprintf('%d ...\n',isub);
    clearvars -except isub roiNames combinedRoiNames
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
    for visualRegion = 1:7
        thisfield = combinedRoiNames{visualRegion};
        load(fullfile(prffolder,['regressPrfSplit' bandpassStr '_v' num2str(visualRegion) '_sub' num2str(isub)  '.mat']), ...
            'nsd', 'numLevels', 'numLengths','rois','nvox','roiPrf','nsplits');
       
        if visualRegion == 4 || visualRegion == 5 || visualRegion == 6 || visualRegion == 7
              
            mean_split1.(thisfield) = squeeze(mean(nsd.voxLenCoef{1}(1,:,1:8),2));
            mean_split2.(thisfield) = squeeze(mean(nsd.voxLenCoef{1}(2,:,1:8),2));

            % average splits    
            nsd.voxLenCoef{1}(nsplits+1,:,:) = mean(nsd.voxLenCoef{1},1);
            mean_all.(thisfield) = squeeze(mean(nsd.voxLenCoef{1}(3,:,1:8),2));
            allLenCoef.(thisfield) = nsd.voxLenCoef{1}(:,:,1:8);

        else %for other regions, combine ventral and dorsal
            % combine ventral and dorsal    
            oldNsd = nsd;
            nsd.voxLenCoef{1} = [];
            nsd.voxLenCoef{2} = [];
            for iroi=1:length(rois)
                nsd.voxLenCoef{1} = cat(2,nsd.voxLenCoef{1},oldNsd.voxLenCoef{iroi});
            end

            mean_split1.(thisfield) = squeeze(mean(nsd.voxLenCoef{1}(1,:,1:8),2));
            mean_split2.(thisfield) = squeeze(mean(nsd.voxLenCoef{1}(2,:,1:8),2));

            % average splits    
            nsd.voxLenCoef{1}(nsplits+1,:,:) = mean(nsd.voxLenCoef{1},1);
            mean_all.(thisfield) = squeeze(mean(nsd.voxLenCoef{1}(3,:,1:8),2));
            allLenCoef.(thisfield) = nsd.voxLenCoef{1}(:,:,1:8);
        end
    end

    saveName = [prffolder, 'voxLenCoef_sub', num2str(isub), '.mat'];
    save(saveName, 'allLenCoef', 'mean_all', 'mean_split1', 'mean_split2');
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
    figure;
    subplot(4,2,1)
    bar(mean_all.V1,1);
    title('v1');

    subplot(4,2,2)
    bar(mean_all.V2,1);
    title('v2');

    subplot(4,2,3)
    bar(mean_all.V3,1);
    title('v3');

    subplot(4,2,4)
    bar(mean_all.hV4,1);
    title('v4');

    subplot(4,2,5)
    bar(mean_all.OPA,1);
    title('OPA');

    subplot(4,2,6)
    bar(mean_all.PPA,1);
    title('PPA');

    subplot(4,2,7)
    bar(mean_all.RSC,1);
    title('RSC');

    totalTitle = ['sub', num2str(isub), ' average for ROI'];
    sgtitle(totalTitle);

    saveas(gcf,['/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Length/len_hist/all_sub', num2str(isub), '.png']);


%% get preferred length for each voxel
% find max coef for each voxel
for visualRegion = 1:7
    thisfield = combinedRoiNames{visualRegion};
    [~, nvox, ~] = size(allLenCoef.(thisfield));
    for ivox = 1: nvox
        [~, ilength] = max(allLenCoef.(thisfield)(3,ivox,:));
        maxCoefLen.(thisfield)(1,ivox) = ilength;
    end
end

% save all ROIs to create overlay
roifolder = ['/bwdata/NSDData/nsddata/ppdata/subj0' num2str(isub) '/func1pt8mm/'];
visualRoisFile = fullfile(roifolder,'roi/prf-visualrois.nii.gz');%V1v, V1d, V2v, V2d, V3v, V3d, and hV4
visRoiData = niftiread(visualRoisFile);
placesRoisFile = fullfile(roifolder,'roi/floc-places.nii.gz'); %OPA, PPA, RSC 
placeRoiData = niftiread(placesRoisFile);

allRoiData = visRoiData;
allRoiData(placeRoiData == 1) = 8;
allRoiData(placeRoiData == 2) = 9;
allRoiData(placeRoiData == 3) = 10;

ourBrain = allRoiData;
ourBrain(ourBrain == 2) = 1;
ourBrain(ourBrain == 3 | ourBrain == 4) = 2;
ourBrain(ourBrain == 5 | ourBrain == 6) = 3;
ourBrain(ourBrain == 7) = 4;
ourBrain(ourBrain == 8) = 5;
ourBrain(ourBrain == 9) = 6;
ourBrain(ourBrain == 10) = 7;

% make a brain volume
newBrain = ourBrain;
newBrain(newBrain > 0) = 0;
for visualRegion = 1:7
    thisfield = combinedRoiNames{visualRegion};
    newBrain(ourBrain == visualRegion) = maxCoefLen.(thisfield)(1,:);
end
newBrain(newBrain < 0) = -1;

for visualRegion = 1:7
    curOurBrain = ourBrain;
    if visualRegion == 3
        curOurBrain(visRoiData == 5 | visRoiData == 6) = 3;
    end
    curNewBrain = curOurBrain;
    curNewBrain(curOurBrain ~= visualRegion) = -1;
    thisfield = combinedRoiNames{visualRegion};

    curNewBrain(curOurBrain == visualRegion) = maxCoefLen.(thisfield)(1,:);

    newBrainbyROI(:,:,:,visualRegion) =curNewBrain;
end

saveName = ['/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Length/brainVolume/lengthBrain_sub', num2str(isub), '.mat'];
save(saveName, 'newBrain');

saveName = ['/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Length/brainVolume/lengthBrainbyROI_sub', num2str(isub), '.mat'];
save(saveName, 'newBrainbyROI');

%% save afni file
% size(newBrain)
cd('/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Length/brainVolume');
%3dinfo betas_session01.nii.gz 
%3dcalc -a betas_session01.nii.gz[1] -expr a -prefix oneBeta
% command = ['3dinfo betas_session01_sub', num2str(isub), '.nii.gz'];
% system(command);
% command = ['3dcalc -a betas_session01_sub', num2str(isub), '.nii.gz[1] -expr a -prefix oneBeta_sub', num2str(isub)];
% system(command);

currOneBeta = ['oneBeta_sub', num2str(isub), '+orig'];
[err,V,Info] = BrikLoad(currOneBeta);

Info.RootName = ['lengthBrain_sub', num2str(isub), '+orig'];
opt.Prefix = ['lengthBrain_sub', num2str(isub)];
WriteBrik(newBrain,Info,opt);
Info.RootName = ['lengthBrainbyROI_sub', num2str(isub), '+orig'];
opt.Prefix = ['lengthBrainbyROI_sub', num2str(isub)];
WriteBrik(newBrainbyROI,Info,opt);

%% save nifti file
% niftiwrite(newBrain,'lengthBrain_sub1.nii');
% info_old = niftiinfo('betas_session01_sub1.nii.gz');
% info_new = niftiinfo('lengthBrain_sub1.nii');
% 
% info_new.PixelDimensions = [1.8, 1.8, 1.8];
% info_new.TransformName = info_old.TransformName;
% info_new.SpatialDimension = info_old.SpatialDimension;
% info_new.Transform = info_old.Transform;
% info_new.Qfactor = info_old.Qfactor;
% info_new.AuxiliaryFile = info_old.AuxiliaryFile;
% info_new.raw.pixdim = info_old.raw.pixdim;
% info_new.raw.aux_file = info_old.raw.aux_file;
% info_new.raw.sform_code = info_old.raw.sform_code;
% info_new.raw.srow_x = info_old.raw.srow_x;
% info_new.raw.srow_y = info_old.raw.srow_y;
% info_new.raw.srow_z = info_old.raw.srow_z;
% 
% niftiwrite(newBrain,'lengthBrain_sub1.nii',info_new);

%% plot figure
figure;
subplot(4,2,1)
histogram(maxCoefLen.V1);
title('v1');

subplot(4,2,2)
histogram(maxCoefLen.V2);
title('v2');

subplot(4,2,3)
histogram(maxCoefLen.V3);
title('v3');

subplot(4,2,4)
histogram(maxCoefLen.hV4);
title('v4');

subplot(4,2,5)
histogram(maxCoefLen.OPA);
title('OPA');

subplot(4,2,6)
histogram(maxCoefLen.PPA);
title('PPA');

subplot(4,2,7)
histogram(maxCoefLen.RSC);
title('RSC');


totalTitle = ['sub', num2str(isub), ' number of preferred length bin'];
sgtitle(totalTitle);

saveas(gcf,['/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Length/len_hist/prefLength_sub', num2str(isub), '.png']);

end
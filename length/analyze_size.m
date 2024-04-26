%analyze_size.m
clear all;
%   uses files created by: analyze_length.m
%   creates files used by:
prffolder = '/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Length/prfsample_Len/';
roiNames = {'V1v','V1d','V2v','V2d','V3v','V3d','hV4','OPA','PPA','RSC'};
combinedRoiNames = {'V1','V2','V3','hV4','OPA','PPA','RSC'};

for isub = 1:8
    clearvars -except isub roiNames combinedRoiNames prffolder

    saveFolder = '/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Length/brainVolume/';
    load([saveFolder, 'lengthBrain_sub', num2str(isub), '.mat']);
    load([saveFolder, 'lengthBrainbyROI_sub', num2str(isub), '.mat']);

    betasfolder = ['/bwdata/NSDData/nsddata/ppdata/subj0' num2str(isub) '/func1pt8mm/'];
    sizeFile = fullfile(betasfolder,'prf_size.nii.gz');
    r2file = fullfile(betasfolder,'prf_R2.nii.gz');
    sizeData = niftiread(sizeFile);
    r2Data = niftiread(r2file);
    sizeData(r2Data<=0) = NaN;
    % image size 8.4° of visual angle. 
    
    % sizeData(sizeData > 50) = NaN;
    sizeData(sizeData > 8.4) = 9;
    sizeData(sizeData >= 7 & sizeData <= 8.4) = 8;
    sizeData(sizeData >= 6 & sizeData < 7) = 7;
    sizeData(sizeData >= 5 & sizeData < 6) = 6;
    sizeData(sizeData >= 4 & sizeData < 5) = 5;
    sizeData(sizeData >= 3 & sizeData < 4) = 4;
    sizeData(sizeData >= 2 & sizeData < 3) = 3;
    sizeData(sizeData >= 1 & sizeData < 2) = 2;
    sizeData(sizeData >= 0 & sizeData < 1) = 1;


    brainsize = size(newBrain);
    sizeBrain = zeros(brainsize(1), brainsize(2),brainsize(3));
    for visualRegion = 1:7
        curNewBrain = newBrainbyROI(:,:,:,visualRegion);
        sizeBrain(curNewBrain ~= -1) = sizeData(curNewBrain ~= -1);
    end
    sizeBrain(sizeBrain == 0) = -1;

    for visualRegion = 1:7
        curNewBrain = newBrainbyROI(:,:,:,visualRegion);
        curNewBrain(curNewBrain~=-1) = sizeData(curNewBrain~=-1);
        sizeBrainbyROI(:,:,:,visualRegion) =curNewBrain;
    end



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


    Info.RootName = ['sizeBrain_sub', num2str(isub), '+orig'];
    opt.Prefix = ['sizeBrain2_sub', num2str(isub)];
    WriteBrik(sizeBrain,Info,opt);


end
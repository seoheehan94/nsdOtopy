% Exclusion: Ecc - 1.177*size >4.2
clear all;
%   uses files created by: analyze_length.m
%   creates files used by:
prffolder = '/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Length/prfsample_Len/';
roiNames = {'V1v','V1d','V2v','V2d','V3v','V3d','hV4','OPA','PPA','RSC'};
combinedRoiNames = {'V1','V2','V3','hV4','OPA','PPA','RSC'};


for isub = 2:8
        clearvars -except isub roiNames combinedRoiNames prffolder

    %% apply exclustion criteria
    
    betasfolder = ['/bwdata/NSDData/nsddata/ppdata/subj0' num2str(isub) '/func1pt8mm/'];
    r2file = fullfile(betasfolder,'prf_R2.nii.gz');
    r2Data = niftiread(r2file);
    eccFile = fullfile(betasfolder,'prf_eccentricity.nii.gz');
    eccData = niftiread(eccFile);
    eccData(r2Data<=0) = NaN;
    sizeFile = fullfile(betasfolder,'prf_size.nii.gz');
    sizeData = niftiread(sizeFile);
    sizeData(r2Data<=0) = NaN;

    exclusionData = eccData - 1.177*sizeData;
    exclusionData(exclusionData>4.2) = NaN;
    % nnz(~isnan(exclusionData))
    % nnz(~isnan(eccData))
    % nnz(~isnan(sizeData))


    %% get median for each ROI
    saveFolder = '/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Length/brainVolume/';
    load([saveFolder, 'lengthBrainbyROI_sub', num2str(isub), '.mat']);

    eccDataNew=eccData;
    eccDataNew(isnan(exclusionData)) = NaN;
    % nnz(~isnan(eccDatanew))

    for visualRegion = 1:7
        curBrain = newBrainbyROI(:,:,:,visualRegion);
        curBrain(curBrain ~= -1) = eccDataNew(curBrain ~= -1);
        curBrain(curBrain == -1) = NaN;
        % eccBrainbyROI(:,:,:,visualRegion) =curNewBrain;

        A = curBrain(~isnan(curBrain));
        topM = curBrain;
        bottomM = curBrain;
        topM(topM<median(A)) = NaN;
        bottomM(bottomM>=median(A)) = NaN;
        
        eccMedBrainbyROI(:,:,:,(2*visualRegion-1)) =topM;
        eccMedBrainbyROI(:,:,:,2*visualRegion) =bottomM;
       

        newBrain = newBrainbyROI(:,:,:,visualRegion);
        topNewBrain = newBrain;
        bottomNewBrain = newBrain;
        topNewBrain(isnan(topM)) = NaN;
        bottomNewBrain(isnan(bottomM)) = NaN;

        newBrainbyROIbyEccMed(:,:,:,(2*visualRegion-1)) =topNewBrain;
        newBrainbyROIbyEccMed(:,:,:,2*visualRegion) =bottomNewBrain;
   
    
    end

    saveName = ['/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Length/brainVolume/eccMedBrainbyROI_sub', num2str(isub), '.mat'];
    % save(saveName, 'newBrainbyROI');

    
   
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


    Info.RootName = ['eccMedBrainbyROI_sub', num2str(isub), '+orig'];
    opt.Prefix = ['eccMedBrainbyROI_sub', num2str(isub)];
    WriteBrik(eccMedBrainbyROI,Info,opt);

    Info.RootName = ['lengthBrainbyROIbyECCMed_sub', num2str(isub), '+orig'];
    opt.Prefix = ['lengthBrainbyROIbyECCMed_sub', num2str(isub)];
    WriteBrik(newBrainbyROIbyEccMed,Info,opt);





end
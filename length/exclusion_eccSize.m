% Exclusion: Ecc - 1.177*size >4.2
clear all;
%   uses files created by: analyze_length.m
%   creates files used by:
prffolder = '/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Length/prfsample_Len/';
roiNames = {'V1v','V1d','V2v','V2d','V3v','V3d','hV4','OPA','PPA','RSC'};
combinedRoiNames = {'V1','V2','V3','hV4','OPA','PPA','RSC'};


for isub = 1:8
        clearvars -except isub roiNames combinedRoiNames prffolder

    %% prf_eccentricity
    saveFolder = '/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Length/brainVolume/';
    load([saveFolder, 'lengthBrain_sub', num2str(isub), '.mat']);
    load([saveFolder, 'lengthBrainbyROI_sub', num2str(isub), '.mat']);

    betasfolder = ['/bwdata/NSDData/nsddata/ppdata/subj0' num2str(isub) '/func1pt8mm/'];
    r2file = fullfile(betasfolder,'prf_R2.nii.gz');
    r2Data = niftiread(r2file);
    eccFile = fullfile(betasfolder,'prf_eccentricity.nii.gz');
    eccData = niftiread(eccFile);
    eccData(r2Data<=0) = NaN;
    sizeFile = fullfile(betasfolder,'prf_size.nii.gz');
    sizeData = niftiread(sizeFile);
    sizeData(r2Data<=0) = NaN;



















end
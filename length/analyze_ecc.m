%analyze_ecc.m
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
    eccFile = fullfile(betasfolder,'prf_eccentricity.nii.gz');
    r2file = fullfile(betasfolder,'prf_R2.nii.gz');
    eccData = niftiread(eccFile);
    r2Data = niftiread(r2file);
    eccData(r2Data<=0) = NaN;
    
    % image size 8.4° of visual angle. Cut off should be half of the size
    % eccData(eccData>4.2) = NaN;
    
    % edges = [0 1 2 3 4.2];
    eccData(eccData > 4.2) = 5;
    eccData(eccData >= 3 & eccData <= 4.2) = 4;
    eccData(eccData >= 2 & eccData < 3) = 3;
    eccData(eccData >= 1 & eccData < 2) = 2;
    eccData(eccData >= 0 & eccData < 1) = 1;

    % get brainvolume
    for curEcc = 1:5
        curNewBrain = newBrain;
        curNewBrain(eccData ~= curEcc) = -1;
        newBrainbyECC(:,:,:,curEcc) =curNewBrain;
    end
 
    for visualRegion = 1:7
        for curEcc = 1:5
            curNewBrain = newBrainbyROI(:,:,:,visualRegion);
            curNewBrain(eccData ~= curEcc) = -1;
            curSubBrik = curEcc + 5*(visualRegion-1);
            newBrainbyROIbyECC(:,:,:,curSubBrik) =curNewBrain;
        end
    end

    brainsize = size(newBrain);
    eccBrain = zeros(brainsize(1), brainsize(2),brainsize(3));
    for visualRegion = 1:7
        curNewBrain = newBrainbyROI(:,:,:,visualRegion);
        eccBrain(curNewBrain ~= -1) = eccData(curNewBrain ~= -1);
    end
    eccBrain(eccBrain == 0) = -1;

    % %% prf-eccrois
    % 
    % eccFile = fullfile(betasfolder,'roi/prf-eccrois.nii.gz');
    % % 0 Unknown
    % % 1 ecc0pt5
    % % 2 ecc1
    % % 3 ecc2
    % % 4 ecc4
    % % 5 ecc4+
    % r2file = fullfile(betasfolder,'prf_R2.nii.gz');
    % eccData = niftiread(eccFile);
    % r2Data = niftiread(r2file);
    % 
    % eccData(eccData == 0 | eccData == -1) = NaN;
    % 
    % % get brainvolume
    % for curEcc = 1:5
    %     curNewBrain = newBrain;
    %     curNewBrain(eccData ~= curEcc) = -1;
    %     newBrainbyECC(:,:,:,curEcc) =curNewBrain;
    % end
    % 
    % for visualRegion = 1:7
    %     for curEcc = 1:5
    %         curNewBrain = newBrainbyROI(:,:,:,visualRegion);
    %         curNewBrain(eccData ~= curEcc) = -1;
    %         curSubBrik = curEcc + 5*(visualRegion-1);
    %         newBrainbyROIbyECC(:,:,:,curSubBrik) =curNewBrain;
    %     end
    % end

    
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


    Info.RootName = ['lengthBrainbyECC_sub', num2str(isub), '+orig'];
    opt.Prefix = ['lengthBrainbyECC_sub', num2str(isub)];
    WriteBrik(newBrainbyECC,Info,opt);

    Info.RootName = ['lengthBrainbyROIbyECC_sub', num2str(isub), '+orig'];
    opt.Prefix = ['lengthBrainbyROIbyECC_sub', num2str(isub)];
    WriteBrik(newBrainbyROIbyECC,Info,opt);


    Info.RootName = ['eccBrain_sub', num2str(isub), '+orig'];
    opt.Prefix = ['eccBrain_sub', num2str(isub)];
    WriteBrik(eccBrain,Info,opt);


   
end
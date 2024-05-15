% Exclusion: Ecc - 1.177*size >4.2
clear all;
%   uses files created by: analyze_length.m
%   creates files used by:
prffolder = '/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Length/prfsample_Len/';
roiNames = {'V1v','V1d','V2v','V2d','V3v','V3d','hV4','OPA','PPA','RSC'};
combinedRoiNames = {'V1','V2','V3','hV4','OPA','PPA','RSC'};

eccDegree = 2;
% allsubMedianValues_ecc=[];
% allsubMedianValues_size=[];

for isub = 1:8
    isub
        clearvars -except isub roiNames combinedRoiNames prffolder allsubMedianValues_ecc allsubMedianValues_size eccDegree

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

    eccDataNew=eccData;
    eccDataNew(isnan(exclusionData)) = NaN;
    sizeDataNew=sizeData;
    sizeDataNew(isnan(exclusionData)) = NaN;
    % nnz(~isnan(eccDatanew))

    % exclude over 30
    eccDataNew(eccDataNew > 30) = NaN;
    sizeDataNew(sizeDataNew > 30) = NaN;

    % AA= unique(eccDataNew(:));
    % AA = AA(~isnan(AA));
    % BB= unique(sizeDataNew(:));
    % BB = BB(~isnan(BB));

    %% load file
    saveFolder = '/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Length/brainVolume/';
    load([saveFolder, 'lengthBrainbyROI_sub', num2str(isub), '.mat']);
    
    %% get median for each ROI
    % medianValues_ecc = zeros(1, 7);
    %  for visualRegion = 1:7
    %     curBrain = newBrainbyROI(:,:,:,visualRegion);
    %     curBrain(curBrain ~= -1) = eccDataNew(curBrain ~= -1);
    %     curBrain(curBrain == -1) = NaN;
    % 
    %     A = curBrain(~isnan(curBrain));
    %     medianValues_ecc(visualRegion) = median(A);
    %     topM = curBrain;
    %     bottomM = curBrain;
    %     topM(topM<median(A)) = NaN;
    %     bottomM(bottomM>=median(A)) = NaN;
    % 
    %     eccMedBrainbyROI(:,:,:,(2*visualRegion-1)) =bottomM;
    %     eccMedBrainbyROI(:,:,:,2*visualRegion) =topM;
    % 
    % 
    %     newBrain = newBrainbyROI(:,:,:,visualRegion);
    %     topNewBrain = newBrain;
    %     bottomNewBrain = newBrain;
    %     topNewBrain(isnan(topM)) = NaN;
    %     bottomNewBrain(isnan(bottomM)) = NaN;
    % 
    %     newBrainbyROIbyEccMed(:,:,:,(2*visualRegion-1)) =bottomNewBrain;
    %     newBrainbyROIbyEccMed(:,:,:,2*visualRegion) =topNewBrain;
    % end
    % 
    % % saveName = ['/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Length/brainVolume/eccMedBrainbyROI_sub', num2str(isub), '.mat'];
    % % save(saveName, 'eccMedBrainbyROI');
    % % 
    % % saveName = ['/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Length/brainVolume/lengthBrainbyROIbyECCMed_sub', num2str(isub), '.mat'];
    % % save(saveName, 'newBrainbyROIbyEccMed');
    % 
    % medianValues_size = zeros(1, 7);
    % for visualRegion = 1:7
    %     curBrain = newBrainbyROI(:,:,:,visualRegion);
    %     curBrain(curBrain ~= -1) = sizeDataNew(curBrain ~= -1);
    %     curBrain(curBrain == -1) = NaN;
    % 
    %     A = curBrain(~isnan(curBrain));
    %     medianValues_size(visualRegion) = median(A);
    %     topM = curBrain;
    %     bottomM = curBrain;
    %     topM(topM<median(A)) = NaN;
    %     bottomM(bottomM>=median(A)) = NaN;
    % 
    %     sizeMedBrainbyROI(:,:,:,(2*visualRegion-1)) =bottomM;
    %     sizeMedBrainbyROI(:,:,:,2*visualRegion) =topM;
    % 
    % 
    %     newBrain = newBrainbyROI(:,:,:,visualRegion);
    %     topNewBrain = newBrain;
    %     bottomNewBrain = newBrain;
    %     topNewBrain(isnan(topM)) = NaN;
    %     bottomNewBrain(isnan(bottomM)) = NaN;
    % 
    %     newBrainbyROIbySizeMed(:,:,:,(2*visualRegion-1)) =bottomNewBrain;
    %     newBrainbyROIbySizeMed(:,:,:,2*visualRegion) =topNewBrain;
    % end
    % 
    % 
    % allsubMedianValues_ecc(isub,:)=medianValues_ecc;
    % % saveName = ['/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Length/brainVolume/sizeMedBrainbyROI_sub', num2str(isub), '.mat'];
    % % save(saveName, 'sizeMedBrainbyROI');
    % % 
    % % saveName = ['/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Length/brainVolume/lengthBrainbyROIbySizeMed_sub', num2str(isub), '.mat'];
    % % save(saveName, 'newBrainbyROIbySizeMed');

    %% center vs. periphery
    
    for visualRegion = 1:7
        curBrain = newBrainbyROI(:,:,:,visualRegion);
        curBrain(curBrain ~= -1) = eccDataNew(curBrain ~= -1);
        curBrain(curBrain == -1) = NaN;

        periph = curBrain;
        center = curBrain;
        periph(periph<eccDegree) = NaN;
        center(center>=eccDegree) = NaN;

        cVSpBrainbyROI(:,:,:,(2*visualRegion-1)) =center;
        cVSpBrainbyROI(:,:,:,2*visualRegion) =periph;


        newBrain = newBrainbyROI(:,:,:,visualRegion);
        periphNewBrain = newBrain;
        centerNewBrain = newBrain;
        periphNewBrain(isnan(periph)) = NaN;
        centerNewBrain(isnan(center)) = NaN;

        newBrainbyROIbycVSp(:,:,:,(2*visualRegion-1)) =centerNewBrain;
        newBrainbyROIbycVSp(:,:,:,2*visualRegion) =periphNewBrain;
    end

    saveName = ['/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Length/brainVolume/eccBrainbyROI_cVSp_sub', num2str(isub), '.mat'];
    save(saveName, 'cVSpBrainbyROI');

    saveName = ['/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Length/brainVolume/lengthBrainbyROIbycVSp_sub', num2str(isub), '.mat'];
    save(saveName, 'newBrainbyROIbycVSp');

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

    Info.RootName = ['eccBrainbyROI_cVSp_sub', num2str(isub), '+orig'];
    opt.Prefix = ['eccBrainbyROI_cVSp_sub', num2str(isub)];
    WriteBrik(cVSpBrainbyROI,Info,opt);

    Info.RootName = ['lengthBrainbyROIbycVSp_sub', num2str(isub), '+orig'];
    opt.Prefix = ['lengthBrainbyROIbycVSp_sub', num2str(isub)];
    WriteBrik(newBrainbyROIbycVSp,Info,opt);

    % 
    % Info.RootName = ['eccMedBrainbyROI_sub', num2str(isub), '+orig'];
    % opt.Prefix = ['eccMedBrainbyROI_sub', num2str(isub)];
    % WriteBrik(eccMedBrainbyROI,Info,opt);
    % 
    % Info.RootName = ['lengthBrainbyROIbyECCMed_sub', num2str(isub), '+orig'];
    % opt.Prefix = ['lengthBrainbyROIbyECCMed_sub', num2str(isub)];
    % WriteBrik(newBrainbyROIbyEccMed,Info,opt);

    % Info.RootName = ['sizeMedBrainbyROI_sub', num2str(isub), '+orig'];
    % opt.Prefix = ['sizeMedBrainbyROI_sub', num2str(isub)];
    % WriteBrik(sizeMedBrainbyROI,Info,opt);
    % 
    % Info.RootName = ['lengthBrainbyROIbySizeMed_sub', num2str(isub), '+orig'];
    % opt.Prefix = ['lengthBrainbyROIbySizeMed_sub', num2str(isub)];
    % WriteBrik(newBrainbyROIbySizeMed,Info,opt);
    % 
    % 



end


%writematrix(allsubMedianValues_ecc,'/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Length/len_hist/allsubMedianValues_ecc.txt');
%writematrix(medianValues_size,'/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Length/len_hist/medianValues_size.txt');

clear all;

isub=1;
visualRegions = [1, 2, 3, 4];
conditions = {'old', 'ori', 'control'};
nsessionsSub = [40 40 32 30 40 32 40 30];
nsessions=nsessionsSub(isub);
nsplits=2;
bandpass = 1; bandMin = 1; bandMax = 7;

boxfolder1 = '/bwdata/NSDData/Seohee/Orientation/prfsample/';
boxfolder2 = '/bwdata/NSDData/Seohee/Orientation/prfsample_Ori/';
savefolder = '/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Orientation/analyses/partition/';
betasfolder = ['/bwdata/NSDData/nsddata_betas/ppdata/subj0' num2str(isub) '/func1pt8mm/betas_fithrf_GLMdenoise_RR/'];
% stimfilename = fullfile(folder,'nsdsynthetic_colorstimuli_subj01.hdf5');
nsdfolder = '/bwdata/NSDData/nsddata/experiments/nsd/';
visualRoisfolder = ['/bwdata/NSDData/nsddata/ppdata/subj0' num2str(isub) '/func1pt8mm/'];

visualRoisFile = fullfile(visualRoisfolder,'roi/prf-visualrois.nii.gz');%V1v, V1d, V2v, V2d, V3v, V3d, and hV4
visRoiData = niftiread(visualRoisFile);
roiNames = {'V1v','V1d','V2v','V2d','V3v','V3d','hV4'};
visRoiData = visRoiData(:);


for visualRegion=visualRegions

    load(fullfile(boxfolder1,['prfSampleStim_v' num2str(visualRegion) '_sub' num2str(isub) '.mat']),'prfSampleLevOri',...
        'rois','allImgs','numLevels','numOrientations','interpImgSize','backgroundSize','pixPerDeg',...
        'roiPrf');
    numLevels1 = numLevels;
    prfSampleLevOri1 = prfSampleLevOri;

    load(fullfile(boxfolder2,['prfSampleStim_ori_v' num2str(visualRegion) '_sub' num2str(isub) '.mat']),'prfSampleLevOri',...
        'rois','allImgs','numLevels','numOrientations','interpImgSize','backgroundSize','pixPerDeg',...
        'roiPrf');
    numLevels2 = numLevels;
    prfSampleLevOri2 = prfSampleLevOri;

    for roinum=1:length(rois); iroi = rois(roinum); roiBetas{roinum}=[]; end
    for isession=1:nsessions
        betasfilename = fullfile(betasfolder,['betas_session' num2str(isession,'%02.f') '.nii.gz']);
        betas = niftiread(betasfilename);
        betas = cast(betas,'double');
        betas = betas/300;
        betas=reshape(betas,[],size(betas,4));
        for roinum=1:length(rois)
            iroi = rois(roinum);
            roiBetas{roinum} = [roiBetas{roinum} betas(visRoiData==iroi,:)];
            roiInd{roinum} = find(visRoiData==iroi);
        end
    end

    nsdDesignFilename = fullfile(nsdfolder, 'nsd_expdesign.mat');
    nsdDesign = load(nsdDesignFilename);
    nsdDesign.masterordering;%for each of 30000 trials, what is corresponding image (out of 10000 images)
    subDesign = nsdDesign.subjectim(isub,nsdDesign.masterordering);%for each of 30000 trials, what is corresponding image (out of 73000 images)
    [imgTrials, imgNum] = ismember(subDesign, allImgs);%logical array

    %if less than 40 sessions, only use image trials that were actually presented
    imgTrials = imgTrials(1:size(roiBetas{roinum},2));
    imgNum = imgNum(1:size(roiBetas{roinum},2));

    imgTrialSum = cumsum(imgTrials);
    splitImgTrials = repmat(imgTrials,2,1);

    midImg = ceil(median(imgNum));
    splitImgTrials(1,imgNum<midImg) = zeros;
    splitImgTrials(2,imgNum>=midImg) = zeros;
    maxNumTrials = max(sum(splitImgTrials,2));

    r2oriSplit1 = cell(length(rois),1);
    r2oriSplit2 = cell(length(rois),1);
    r2oriSplitCombined = cell(length(rois),1);



    for roinum=1:length(rois)
        nvox(roinum) = size(roiBetas{roinum},1);
        voxOriResidualSplit1{roinum} = NaN(nsplits, nvox(roinum),maxNumTrials);
        voxOriResidualSplit2{roinum} = NaN(nsplits, nvox(roinum),maxNumTrials);
        voxOriResidualSplitCombined{roinum} = NaN(nsplits, nvox(roinum),maxNumTrials);
        voxOriCoef1{roinum} = zeros(nsplits, nvox(roinum),numLevels1*numOrientations+1);
        voxOriCoef2{roinum} = zeros(nsplits, nvox(roinum),numLevels2*numOrientations+1);
        voxOriCoefCombined{roinum} = zeros(nsplits, nvox(roinum),numLevels1*numOrientations+numLevels2*numOrientations+1);
        
        for isplit=1:nsplits
            imgTrials = splitImgTrials(isplit,:);
            numTrials = sum(imgTrials);
            for ivox=1:nvox(roinum)
                voxBetas = roiBetas{roinum}(ivox,imgTrials>0)';

                voxPrfOriSample1 = squeeze(prfSampleLevOri1{roinum}(imgNum(imgTrials>0),ivox,:,:));
                voxPrfOriSample1 = reshape(voxPrfOriSample1,[],numLevels1*numOrientations);
                voxPrfOriSample2 = squeeze(prfSampleLevOri2{roinum}(imgNum(imgTrials>0),ivox,:,:));
                voxPrfOriSample2 = reshape(voxPrfOriSample2,[],numLevels2*numOrientations);
                voxPrfOriSampleCombined = [voxPrfOriSample1, voxPrfOriSample2];

                %add constant predictor
                voxPrfOriSample1(:,end+1) = ones;
                voxPrfOriSample2(:,end+1) = ones;
                voxPrfOriSampleCombined(:,end+1) = ones;

                voxOriCoef1{roinum}(isplit,ivox,:) = voxPrfOriSample1\voxBetas;%check vox 144 in first ROI
                voxOriCoef2{roinum}(isplit,ivox,:) = voxPrfOriSample2\voxBetas;
                voxOriCoefCombined{roinum}(isplit,ivox,:) = voxPrfOriSampleCombined\voxBetas;
            end
        end

        for isplit=1:nsplits
            imgTrials = splitImgTrials(isplit,:);
            numTrials = sum(imgTrials);
            for ivox=1:nvox(roinum)

                voxBetas = roiBetas{roinum}(ivox,imgTrials>0)';

                voxPrfOriSample1 = squeeze(prfSampleLevOri1{roinum}(imgNum(imgTrials>0),ivox,:,:));
                voxPrfOriSample1 = reshape(voxPrfOriSample1,[],numLevels1*numOrientations);
                voxPrfOriSample2 = squeeze(prfSampleLevOri2{roinum}(imgNum(imgTrials>0),ivox,:,:));
                voxPrfOriSample2 = reshape(voxPrfOriSample2,[],numLevels2*numOrientations);
                voxPrfOriSampleCombined = [voxPrfOriSample1, voxPrfOriSample2];

                %add constant predictor
                voxPrfOriSample1(:,end+1) = ones;
                voxPrfOriSample2(:,end+1) = ones;
                voxPrfOriSampleCombined(:,end+1) = ones;

                voxOriResidualSplit1{roinum}(isplit,ivox,1:numTrials) = voxBetas' - squeeze(voxOriCoef1{roinum}(nsplits-isplit+1,ivox,:))'*voxPrfOriSample1';
                voxOriResidualSplit2{roinum}(isplit,ivox,1:numTrials) = voxBetas' - squeeze(voxOriCoef2{roinum}(nsplits-isplit+1,ivox,:))'*voxPrfOriSample2';
                voxOriResidualSplitCombined{roinum}(isplit,ivox,1:numTrials) = voxBetas' - squeeze(voxOriCoefCombined{roinum}(nsplits-isplit+1,ivox,:))'*voxPrfOriSampleCombined';

                %r2 between splits
                r2oriSplit1{roinum}(isplit,ivox) = rsquared(voxOriResidualSplit1{roinum}(isplit,ivox,1:sum(splitImgTrials(isplit,:))), roiBetas{roinum}(ivox,imgTrials>0));
                r2oriSplit2{roinum}(isplit,ivox) = rsquared(voxOriResidualSplit2{roinum}(isplit,ivox,1:sum(splitImgTrials(isplit,:))), roiBetas{roinum}(ivox,imgTrials>0));
                r2oriSplitCombined{roinum}(isplit,ivox) = rsquared(voxOriResidualSplitCombined{roinum}(isplit,ivox,1:sum(splitImgTrials(isplit,:))), roiBetas{roinum}(ivox,imgTrials>0));

                % Partition variance
                sharedVariance{roinum}(isplit,ivox) = r2oriSplitCombined{roinum}(isplit,ivox) - (r2oriSplit1{roinum}(isplit,ivox) + r2oriSplit2{roinum}(isplit,ivox));
                uniqueModel1{roinum}(isplit,ivox) = r2oriSplit1{roinum}(isplit,ivox) - sharedVariance{roinum}(isplit,ivox);
                uniqueModel2{roinum}(isplit,ivox) = r2oriSplit2{roinum}(isplit,ivox) - sharedVariance{roinum}(isplit,ivox);
            end
        end
    end

    nsd.r2oriSplitCombined = r2oriSplitCombined;
    nsd.sharedVariance = sharedVariance;
    nsd.uniqueModel1 = uniqueModel1;
    nsd.uniqueModel2 = uniqueModel2;


     save(fullfile(savefolder,['regressPartition_oldori_v' num2str(visualRegion) '_sub' num2str(isub) '.mat']), ...
        'nsd','rois','nvox','nsplits');


end
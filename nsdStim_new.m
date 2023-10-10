% nsdStim.m
%
% associated with the following publication: Roth, ZN, Kay, K, and Merriam, EP (2022).
% Massive natural scene sampling reveals reliable coarse-scale orientation tuning in human V1
% DOI:
%
%   usage: nsdStim()
%   by: zvi roth
%   date: 7/29/2022
%   purpose: run natural scene stimuli through the models, get filter energy responses,
%   and save output
%   creates files used by: prfSampleModel.m

% uses the steerable pyramid: https://github.com/elimerriam/stimulusVignetting


close all
clear all

interpImgSize = 714;
backgroundSize = 1024;
imgScaling = 0.5;

%pyramidfolder = '/misc/data18/rothzn/nsd/stimuli/pyramid/';%to save model outputs
orifolder = '/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/stimuli/orientationfilter/';%to save model outputs

%%
% construct quad frequency filters
numOrientations = 8;
bandwidth = 1;
dims = [backgroundSize backgroundSize];
dims = dims*imgScaling;
numLevels = 1;
[freqRespsImag, freqRespsReal, pind] = makeQuadFRs(dims, numLevels, numOrientations, bandwidth);

%%
stimfolder = '/bwdata/NSDData/stmuli/';
stimfilename = fullfile(stimfolder,'nsd_stimuli.hdf5');%[3 1360 714 220]
stiminfo = h5info(stimfilename);

imgSizeX = stiminfo.Datasets.Dataspace.Size(2);%1360;
imgSizeY = stiminfo.Datasets.Dataspace.Size(3);%714;
numImgs = stiminfo.Datasets.Dataspace.Size(4);%220

nsdfolder = '/bwdata/NSDData/nsddata/experiments/nsd/';
nsdDesignFilename = fullfile(nsdfolder, 'nsd_expdesign.mat');
nsdDesign = load(nsdDesignFilename);
allImgs = nsdDesign.sharedix; %indices of the shared 1000 images

vecLDfolder = '/bwlab/Users/SeoheeHan/NSDData/nsddata_stimuli';

for isub=[1:1]
    
    
    allImgs = nsdDesign.subjectim(isub,nsdDesign.masterordering);%indices of all 10000 images used for this subject
    allImgs = unique(allImgs);
    %%
    backgroundColor(1,1,:) = uint8([127,127,127]);
    
    nImgs = 1;
    
    fixPoint(1,1,:) = [255 0 0];
    
    iimg=0;
    for imgNum=allImgs
        iimg = iimg+1
        
        orifilename = ['oriImg' num2str(imgNum) '.mat'];
        if ~isfile(fullfile(orifolder, orifilename))%if file exists already no need to remake it
            origImg = h5read(stimfilename,'/imgBrick/',[1 1 1 imgNum],[3 imgSizeX imgSizeY nImgs]);
            origImg = double(origImg);
            origImg = permute(origImg,[3 2 1]);%[425,425,3]
            
            [Xq, Yq] = meshgrid(linspace(1,imgSizeX, interpImgSize), linspace(1,imgSizeY, interpImgSize));
            for irgb=1:3
                interpImg(:,:,irgb) = interp2(squeeze(origImg(:,:,irgb)), Xq, Yq);
            end
            %add red semi-transparent fixation point
            
            interpImg(interpImgSize/2-8:interpImgSize/2+8,interpImgSize/2-8:interpImgSize/2+8,:) = ...
                (interpImg(interpImgSize/2-8:interpImgSize/2+8,interpImgSize/2-8:interpImgSize/2+8,:) + ...
                repmat(fixPoint,17,17,1))/2;
            
            %%add background
            bigImg = repmat(backgroundColor,backgroundSize,backgroundSize,1);
            bigImg(1+backgroundSize/2-interpImgSize/2 : backgroundSize/2+interpImgSize/2, 1+backgroundSize/2-interpImgSize/2 : backgroundSize/2+interpImgSize/2,:) = interpImg(:,:,:);
            
            %axis image
            % colormap gray
            
            %change to grayscale
            % for now, simply by averaging across RGB channels
            bigImg = mean(bigImg,3);
            
            %DOWNSAMPLE
            bigImg = imresize(bigImg,imgScaling);
            %% pass image through orientation filter
            oriMap = generateOrientationMap(vecLD, NaN, [512 512]);
            binWidth2 = 90 / length(vecLD.orientationBins);
            horIdx = (oriMap > (180-binWidth2));
            oriMap(horIdx) = oriMap(horIdx) - 180;
            
            %vecLD.orientationBins
            for binIdx = 1: length(vecLD.orientationBins)
                oribinmap{1,binIdx} = (abs(oriMap-vecLD.orientationBins(binIdx)) <= binWidth2);

            end
            modelOri=cat(3,oribinmap{:});
            modelOri = permute(modelOri,[3 1 2]);

            save(fullfile(orifolder, orifilename), 'interpImgSize','backgroundSize','imgScaling',...
                'numOrientations','bandwidth','dims','bigImg','sumOri','modelOri','numLevels');
        end
    end
end


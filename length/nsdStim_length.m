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

close all;
clear all;

backgroundSize = [512 512];
renderSize = [357,357];

addpath(genpath('/home/hanseohe/Documents/GitHub/nsdOtopy'));
addpath(genpath('/home/hanseohe/Documents/GitHub/mrTools'));
addpath(genpath('/home/hanseohe/Documents/GitHub/stimulusVignetting'));

%cd '/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/stimuli';
%pyramidfolder = '/misc/data18/rothzn/nsd/stimuli/pyramid/';%to save model outputs
lenfolder = '/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/stimuli/lengthfilter/';%to save model outputs

%%
% construct quad frequency filters

numLengths = 8;
bandwidth = 1;
dims = backgroundSize;
numLevels = 1;
[freqRespsImag, freqRespsReal, pind] = makeQuadFRs(dims, numLevels, numLengths, bandwidth);

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

for isub=[1:8]
    
    allImgs = nsdDesign.subjectim(isub,nsdDesign.masterordering);%indices of all 10000 images used for this subject
    allImgs = unique(allImgs);
    %%
    backgroundColor(1,1,:) = uint8([127,127,127]);
    
    nImgs = 1;
    
    fixPoint(1,1,:) = [255 0 0];
    
    iimg=0;
    for imgNum=allImgs
        iimg = iimg+1
        
        filename = ['lenImg' num2str(imgNum) '.mat'];
        if ~isfile(fullfile(lenfolder, filename))%if file exists already no need to remake it
            imgName = ['img' num2str(imgNum) '.mat'];
            if imgNum <= 14600*1
                imgFolder = '/bwlab/Users/SeoheeHan/NSDData/nsddata_stimuli/images01/';
            elseif imgNum <= 14600*2
                imgFolder = '/bwlab/Users/SeoheeHan/NSDData/nsddata_stimuli/images02/';
            elseif imgNum <= 14600*3
                imgFolder = '/bwlab/Users/SeoheeHan/NSDData/nsddata_stimuli/images03/';
            elseif imgNum <= 14600*4
                imgFolder = '/bwlab/Users/SeoheeHan/NSDData/nsddata_stimuli/images04/';
            elseif imgNum <= 14600*5
                imgFolder = '/bwlab/Users/SeoheeHan/NSDData/nsddata_stimuli/images05/';
            end

            load(fullfile(imgFolder, imgName));
           
            %% pass image through orientation filter
            lenMap = generateLengthMap(vecLD, NaN, backgroundSize, renderSize);

            minmaxLength = [2,sum(vecLD.imsize)];
            logMinMax = log10(minmaxLength + 1);
            numBins =8;
            binWidth = (logMinMax(2)-logMinMax(1)) / numBins; %the range of the original length is from max to min length value
            binBoundary = [logMinMax(1) : binWidth : logMinMax(2)];
            binMax = 10.^binBoundary - 1;
            
            %vecLD.orientationBins
            for binIdx = 1: length(vecLD.lengthBins)
                lenbinmap{1,binIdx} = (binMax(binIdx) <= lenMap) & (lenMap <= binMax(binIdx+1));
            end
            modelLen=cat(3,lenbinmap{:});
            modelLen = permute(modelLen,[3 1 2]);

            save(fullfile(lenfolder, filename),...
                'numLengths','bandwidth','dims','modelLen','numLevels', 'lenMap');
        end
    end
end


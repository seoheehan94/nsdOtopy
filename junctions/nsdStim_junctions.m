% nsdStim_junctions.m
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

cd '/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/stimuli';
%pyramidfolder = '/misc/data18/rothzn/nsd/stimuli/pyramid/';%to save model outputs
juncfolder = '/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/stimuli/juncfilter/';%to save model outputs

%%
% construct quad frequency filters

numJunctions = 4;
bandwidth = 1;
dims = backgroundSize;
numLevels = 1;
[freqRespsImag, freqRespsReal, pind] = makeQuadFRs(dims, numLevels, numJunctions, bandwidth);

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
        
        juncfilename = ['juncImg' num2str(imgNum) '.mat'];
        if ~isfile(fullfile(juncfolder, juncfilename))%if file exists already no need to remake it
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
           
            %% pass image through junction filter
            juncMap = generateJunctionMap(vecLD, NaN, backgroundSize, renderSize);
   
            binWidth2 = 90 / length(vecLD.orientationBins);
            horIdx = (juncMap > (180-binWidth2));
            juncMap(horIdx) = juncMap(horIdx) - 180;
            
            %vecLD.orientationBins
            for binIdx = 1: length(vecLD.orientationBins)
                juncbinmap{1,binIdx} = (abs(juncMap-vecLD.orientationBins(binIdx)) <= binWidth2);

            end
            modelJunc=cat(3,juncbinmap{:});
            modelJunc = permute(modelJunc,[3 1 2]);

            save(fullfile(juncfolder, juncfilename),...
                'numJunctions','bandwidth','dims','modelJunc','numLevels', 'juncMap');
        end
    end
end


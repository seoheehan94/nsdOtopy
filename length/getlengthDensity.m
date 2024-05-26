% getlengthDensity.m
clear all;

imgFolder = '/bwlab/Users/SeoheeHan/NSDData/nsddata_stimuli/';
imgSubFolder=dir(imgFolder); 
imgSubFolder = imgSubFolder(ismember({imgSubFolder(:).name},{'images01','images02','images03','images04','images05'}));

load("/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Length/Centervs.Periphery/centerMask.mat");
load("/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Length/Centervs.Periphery/peripheryMask.mat");

for k=3:length(imgSubFolder)
    thisFolder = imgSubFolder(k).name;
    imFilePath = natsortfiles(dir([imgFolder, thisFolder,'/*.mat']));

    normSumList = [];
    for j = 1:length(imFilePath)
        imFile = imFilePath(j).name;
        fprintf('%d. %d. %s ...\n',k,j,imFile);
        load([imgFolder, thisFolder, '/', imFile]);

        %run generateFeatureDensityMap
        FDM = generateFeatureDensityMap(vecLD,'length');

        %Apply mask center/periphery (about the same number of pixels)
        FDMcenter = FDM;
        FDMcenter(centerMask==0) =0;
        FDMperiph = FDM;
        FDMperiph(peripheryMask==0) =0;

        figure;imagesc(FDMcenter)
        figure;imagesc(FDMperiph)
        %Sum all the values divided by non-zero pixels
        sumFDMcenter = sum(FDMcenter(:));
        numPixcenter = length(FDMcenter(FDMcenter~=0));
        normSumcenter=sumFDMcenter/numPixcenter;

        sumFDMperiph = sum(FDMperiph(:));
        numPixperiph = length(FDMperiph(FDMperiph~=0));
        normSumperiph=sumFDMperiph/numPixperiph;

        normSumList(j,:)=[normSumcenter,normSumperiph];

        
    end

    save(['/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Length/Centervs.Periphery/normFDMSum0',num2str(k)],'normSumList');
 end




%%
fileFolder = '/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Length/Centervs.Periphery/';
fileDir = dir([fileFolder, 'norm*']);
TotalList = [];
for k=1:length(fileDir)
    load([fileFolder, fileDir(k).name]);
    TotalList = [TotalList; normSumList];
end

%center, periphery
mean(TotalList, "omitnan")
median(TotalList, "omitnan")


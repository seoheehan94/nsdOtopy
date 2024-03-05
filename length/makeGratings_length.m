%makeGratings_length.m
clear all;
rng(4228);


%% check average number of pixels of images
% original_folder = '/bwlab/Users/SeoheeHan/NSDData/nsddata_stimuli'; % Original image file folder
% addpath(original_folder);
% original_subfolder = dir(original_folder);
% %original_subfolder = original_subfolder(~ismember({original_subfolder(:).name},{'.','..'}));
% original_subfolder = original_subfolder(ismember({original_subfolder(:).name},{'images01','images02','images03','images04','images05'}));
% 
% meanPixelList = [];
% 
% for k = 1 : length(original_subfolder)
%     name = original_subfolder(k).name;
%     imFilePath = natsortfiles(dir(strcat(original_folder, '/', name,'/', '*.mat')));
%     vecLD = [];
% 
%     for j = 1 : length(imFilePath)
%         imFile = imFilePath(j).name;
%         fprintf('%d. %s ...\n',j,imFile);
%         load([original_folder, '/', name, '/', imFile])
%         meanPixelList = [meanPixelList;sum(vecLD.contourLengths)];
%     end
% end
% 
% mean(meanPixelList)

% grandMean = 6688.66;
grandMean = 6689;

%% get contours from one length bin
for curbin = 8:8
    fprintf('%d\n',curbin);
    % curbin = 1;
    imgOrder = 1:73000;
    imgOrder = imgOrder(randperm(length(imgOrder)));
    newVecLD = struct;
    testvecLD = struct;
    newVecLD.contours =[];
    newVecLD.numContours = 0;
    newVecLD.contourLengths = [];
    imgList = [];
    
    for k = 1: length(imgOrder)
        % load vecLD file
        vecLD = [];
        imgName = ['img' num2str(imgOrder(k)) '.mat'];
        if imgOrder(k) <= 14600*1
            imgFolder = '/bwlab/Users/SeoheeHan/NSDData/nsddata_stimuli/images01/';
        elseif imgOrder(k) <= 14600*2
            imgFolder = '/bwlab/Users/SeoheeHan/NSDData/nsddata_stimuli/images02/';
        elseif imgOrder(k) <= 14600*3
            imgFolder = '/bwlab/Users/SeoheeHan/NSDData/nsddata_stimuli/images03/';
        elseif imgOrder(k) <= 14600*4
            imgFolder = '/bwlab/Users/SeoheeHan/NSDData/nsddata_stimuli/images04/';
        elseif imgOrder(k) <= 14600*5
            imgFolder = '/bwlab/Users/SeoheeHan/NSDData/nsddata_stimuli/images05/';
        end
        load(fullfile(imgFolder, imgName));

        % get contour list for the current histogram bin
        curbinList = vecLD.lengthHistograms(:,curbin);
        curIdxList = find(curbinList);

        % if the list is not empty
        if ~isempty(curIdxList)
            done = 0;
            c = 1;
            while done == 0
                % randomly choose one contour
                randCurIdxList = curIdxList(randperm(length(curIdxList)));

                thisIdx = randCurIdxList(c);
                thisContour = vecLD.contours(thisIdx);

                % test if the new contour intersect with previous contours
                testvecLD = newVecLD;
                testvecLD.contours = [testvecLD.contours, thisContour];
                testvecLD.numContours = testvecLD.numContours+1;
                testvecLD.contourLengths = [testvecLD.contourLengths, vecLD.contourLengths(thisIdx)];

                Junctions = detectJunctions(testvecLD);
                if isempty(Junctions) % if no junctions add the contour
                    newVecLD = testvecLD;
                    imgList = [imgList; {imgName, thisIdx}];

                    done = 1;
                else % if yes junctions move on to the next contour
                    c = c+1;
                    if c > length(curIdxList)
                        done = 1;
                    end
                end
            end
        end

        % if the total number of pixels is larger than the mean, stop
        if sum(newVecLD.contourLengths) > grandMean
            break
        end

    end
    % save newVecLD
    newVecLD.imsize = vecLD.imsize;
    newVecLD.imgList = imgList;
    save(['/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Length/grating/gratingv4_', num2str(curbin), '.mat'], 'newVecLD');

    % save as an image
    imgLD = renderLinedrawing(newVecLD);
    imgLD = squeeze(imgLD(:,:,1)); % use grayscale encoding
    imwrite(imgLD,['/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Length/grating/gratingv4_', num2str(curbin), '.png']);
end

%% pass image through orientation filter
backgroundSize = [512 512];
renderSize = [357,357];

numLengths = 8;
bandwidth = 1;
dims = backgroundSize;
numLevels = 1;
[freqRespsImag, freqRespsReal, pind] = makeQuadFRs(dims, numLevels, numLengths, bandwidth);

imgFolder='/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Length/grating/'; 
for imgNum = 1:8
    imgName = ['grating' num2str(imgNum) '.mat'];
    load(fullfile(imgFolder, imgName));
    newVecLD = computeContourProperties(newVecLD);
    newVecLD = getContourPropertiesStats(newVecLD);

    lenMap = generateLengthMap(newVecLD, NaN, backgroundSize, renderSize);

    minmaxLength = [2,sum(newVecLD.imsize)];
    logMinMax = log10(minmaxLength + 1);
    numBins =8;
    binWidth = (logMinMax(2)-logMinMax(1)) / numBins; %the range of the original length is from max to min length value
    binBoundary = [logMinMax(1) : binWidth : logMinMax(2)];
    logbins = binBoundary(2:end) - binWidth/2;
    %binMax = 10.^binBoundary - 1;

    %vecLD.orientationBins
    for binIdx = 1: length(newVecLD.lengthBins)
        lenbinmap{1,binIdx} = (abs(lenMap-logbins(binIdx)) <= binWidth/2);

        
        modelLenEnergy(1,imgNum,1,binIdx) = sum(lenbinmap{1,binIdx}(:));
    end

    sumLenEnergy(1,imgNum,1) = sum(modelLenEnergy(1,imgNum,1,:));
end



save('gratings_length.mat','numLengths','numLevels',...
    'sumLenEnergy','modelLenEnergy');






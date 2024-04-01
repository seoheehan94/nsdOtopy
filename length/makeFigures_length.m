% makeFigures_length.m
clear all;
cd('/home/hanseohe/Documents/GitHub/nsdOtopy/length');

prffolder = '/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Length/prfsample_Len/';  
roiNames = {'V1v','V1d','V2v','V2d','V3v','V3d','hV4','OPA','PPA','RSC'};
combinedRoiNames = {'V1','V2','V3','hV4','OPA','PPA','RSC'};

%% 1.  Coef voxel preference %%
lenPref = struct;
fullecc_lenpref = struct;
ecc_lenpref = struct;
for isub = 1:8
    % load(fullfile([prffolder, 'voxLenCoef_sub', num2str(isub), '.mat']));
 
    % subName = ['sub', num2str(isub)];
    % voxPref.(subName)=maxCoefLen;
    %% check consistency between splits
    % for visualRegion = 1:4
    %     thisfield = ['v', num2str(visualRegion)];
    %
    %     figure;
    %     subplot(1,2,1)
    %     bar(mean_split1.(thisfield),1);
    %     title('split1');
    %
    %     subplot(1,2,2)
    %     bar(mean_split2.(thisfield),1);
    %     title('split2');
    %
    %     totalTitle = ['sub', num2str(isub), ' V', num2str(visualRegion)];
    %     sgtitle(totalTitle);
    %
    %     saveas(gcf,['/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Length/len_hist/split_sub', num2str(isub), 'v', num2str(visualRegion), '.png']);
    % end

    %% get average for each ROI
    % figure;
    % subplot(4,2,1)
    % bar(mean_all.V1,1);
    % title('v1');
    % 
    % subplot(4,2,2)
    % bar(mean_all.V2,1);
    % title('v2');
    % 
    % subplot(4,2,3)
    % bar(mean_all.V3,1);
    % title('v3');
    % 
    % subplot(4,2,4)
    % bar(mean_all.hV4,1);
    % title('v4');
    % 
    % subplot(4,2,5)
    % bar(mean_all.OPA,1);
    % title('OPA');
    % 
    % subplot(4,2,6)
    % bar(mean_all.PPA,1);
    % title('PPA');
    % 
    % subplot(4,2,7)
    % bar(mean_all.RSC,1);
    % title('RSC');
    % 
    % totalTitle = ['sub', num2str(isub), ' average for ROI'];
    % sgtitle(totalTitle);
    % 
    % saveas(gcf,['/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Length/len_hist/all_sub', num2str(isub), '.png']);

    %% preferred length from coef
    % figure;
    % subplot(4,2,1)
    % histogram(maxCoefLen.V1);
    % title('v1');
    % 
    % subplot(4,2,2)
    % histogram(maxCoefLen.V2);
    % title('v2');
    % 
    % subplot(4,2,3)
    % histogram(maxCoefLen.V3);
    % title('v3');
    % 
    % subplot(4,2,4)
    % histogram(maxCoefLen.hV4);
    % title('v4');
    % 
    % subplot(4,2,5)
    % histogram(maxCoefLen.OPA);
    % title('OPA');
    % 
    % subplot(4,2,6)
    % histogram(maxCoefLen.PPA);
    % title('PPA');
    % 
    % subplot(4,2,7)
    % histogram(maxCoefLen.RSC);
    % title('RSC');
    % 
    % 
    % totalTitle = ['sub', num2str(isub), ' number of preferred length bin'];
    % sgtitle(totalTitle);
    % 
    % saveas(gcf,['/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Length/len_hist/prefLength_sub', num2str(isub), '.png']);

    %% Eccentricity scatterplot
    clearvars eccBrainbyROI 
    
    saveFolder = '/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Length/brainVolume/';
    load([saveFolder, 'lengthBrain_sub', num2str(isub), '.mat']);
    load([saveFolder, 'lengthBrainbyROI_sub', num2str(isub), '.mat']);

    betasfolder = ['/bwdata/NSDData/nsddata/ppdata/subj0' num2str(isub) '/func1pt8mm/'];
    eccFile = fullfile(betasfolder,'prf_eccentricity.nii.gz');
    r2file = fullfile(betasfolder,'prf_R2.nii.gz');
    eccData = niftiread(eccFile);
    r2Data = niftiread(r2file);
    eccData(r2Data<=0) = NaN;

    % eccentricity < 20
    for visualRegion = 1:7
        curNewBrain = newBrainbyROI(:,:,:,visualRegion);
        curNewBrain(curNewBrain~=-1) = eccData(curNewBrain~=-1);
        eccBrainbyROI(:,:,:,visualRegion) =curNewBrain;
    end
    for visualRegion = 1:7
        thisfield = ['sub', num2str(isub)];
        thisEcc = eccBrainbyROI(:,:,:,visualRegion);
        thisEcc = thisEcc(thisEcc ~= -1);
        thisVoxPref = newBrainbyROI(:,:,:,visualRegion);
        thisVoxPref = thisVoxPref(thisVoxPref ~= -1);
        fullecc_lenpref.(thisfield){visualRegion}(:,1)=thisEcc;
        fullecc_lenpref.(thisfield){visualRegion}(:,2)=thisVoxPref;
        [row,col]=find(fullecc_lenpref.(thisfield){visualRegion}>=20);
        fullecc_lenpref.(thisfield){visualRegion}(row,col)= NaN;
    end

    % figure;
    % subplot(4,2,1)
    % scatter(fullecc_lenpref.(thisfield){1}(:,2), fullecc_lenpref.(thisfield){1}(:,1))
    % [r,p] =corr(fullecc_lenpref.(thisfield){1}(:,2), fullecc_lenpref.(thisfield){1}(:,1),'rows','complete');
    % title(['v1: ', num2str(r), '/', num2str(p)]);
    % subplot(4,2,2)
    % scatter(fullecc_lenpref.(thisfield){2}(:,2), fullecc_lenpref.(thisfield){2}(:,1))
    % [r,p] =corr(fullecc_lenpref.(thisfield){2}(:,2), fullecc_lenpref.(thisfield){2}(:,1),'rows','complete');
    % title(['v2: ', num2str(r), '/', num2str(p)]);
    % subplot(4,2,3)
    % scatter(fullecc_lenpref.(thisfield){3}(:,2), fullecc_lenpref.(thisfield){3}(:,1))
    % [r,p] =corr(fullecc_lenpref.(thisfield){3}(:,2), fullecc_lenpref.(thisfield){3}(:,1),'rows','complete');
    % title(['v3: ', num2str(r), '/', num2str(p)]);
    % subplot(4,2,4)
    % scatter(fullecc_lenpref.(thisfield){4}(:,2), fullecc_lenpref.(thisfield){4}(:,1))
    % [r,p] =corr(fullecc_lenpref.(thisfield){4}(:,2), fullecc_lenpref.(thisfield){4}(:,1),'rows','complete');
    % title(['v4: ', num2str(r), '/', num2str(p)]);
    % subplot(4,2,5)
    % scatter(fullecc_lenpref.(thisfield){5}(:,2), fullecc_lenpref.(thisfield){5}(:,1))
    % [r,p] =corr(fullecc_lenpref.(thisfield){5}(:,2), fullecc_lenpref.(thisfield){5}(:,1),'rows','complete');
    % title(['OPA: ', num2str(r), '/', num2str(p)]);
    % subplot(4,2,6)
    % scatter(fullecc_lenpref.(thisfield){6}(:,2), fullecc_lenpref.(thisfield){6}(:,1))
    % [r,p] =corr(fullecc_lenpref.(thisfield){6}(:,2), fullecc_lenpref.(thisfield){6}(:,1),'rows','complete');
    % title(['PPA: ', num2str(r), '/', num2str(p)]);
    % subplot(4,2,7)
    % scatter(fullecc_lenpref.(thisfield){7}(:,2), fullecc_lenpref.(thisfield){7}(:,1))
    % [r,p] =corr(fullecc_lenpref.(thisfield){7}(:,2), fullecc_lenpref.(thisfield){7}(:,1),'rows','complete');
    % title(['RSC: ', num2str(r), '/', num2str(p)]);
    % totalTitle = ['sub', num2str(isub), ' preferred length x full eccentricity'];
    % sgtitle(totalTitle);

    % saveas(gcf,['/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Length/len_hist/prefLen_fullE', num2str(isub), '.png']);

    % eccentricity <= 4.2
    eccData(eccData > 4.2) = NaN;
    for visualRegion = 1:7
        curNewBrain = newBrainbyROI(:,:,:,visualRegion);
        curNewBrain(curNewBrain~=-1) = eccData(curNewBrain~=-1);
        eccBrainbyROI(:,:,:,visualRegion) =curNewBrain;
    end
    for visualRegion = 1:7
        thisfield = ['sub', num2str(isub)];
        thisEcc = eccBrainbyROI(:,:,:,visualRegion);
        thisEcc = thisEcc(thisEcc ~= -1);
        thisVoxPref = newBrainbyROI(:,:,:,visualRegion);
        thisVoxPref = thisVoxPref(thisVoxPref ~= -1);
        ecc_lenpref.(thisfield){visualRegion}(:,1)=thisEcc;
        ecc_lenpref.(thisfield){visualRegion}(:,2)=thisVoxPref;
    end
    
    % figure;
    % subplot(4,2,1)
    % scatter(ecc_lenpref.(thisfield){1}(:,2), ecc_lenpref.(thisfield){1}(:,1))
    % [r,p] =corr(ecc_lenpref.(thisfield){1}(:,2), ecc_lenpref.(thisfield){1}(:,1),'rows','complete');
    % title(['v1: ', num2str(r), '/', num2str(p)]);
    % subplot(4,2,2)
    % scatter(ecc_lenpref.(thisfield){2}(:,2), ecc_lenpref.(thisfield){2}(:,1))
    % [r,p] =corr(ecc_lenpref.(thisfield){2}(:,2), ecc_lenpref.(thisfield){2}(:,1),'rows','complete');
    % title(['v2: ', num2str(r), '/', num2str(p)]);
    % subplot(4,2,3)
    % scatter(ecc_lenpref.(thisfield){3}(:,2), ecc_lenpref.(thisfield){3}(:,1))
    % [r,p] =corr(ecc_lenpref.(thisfield){3}(:,2), ecc_lenpref.(thisfield){3}(:,1),'rows','complete');
    % title(['v3: ', num2str(r), '/', num2str(p)]);
    % subplot(4,2,4)
    % scatter(ecc_lenpref.(thisfield){4}(:,2), ecc_lenpref.(thisfield){4}(:,1))
    % [r,p] =corr(ecc_lenpref.(thisfield){4}(:,2), ecc_lenpref.(thisfield){4}(:,1),'rows','complete');
    % title(['v4: ', num2str(r), '/', num2str(p)]);
    % subplot(4,2,5)
    % scatter(ecc_lenpref.(thisfield){5}(:,2), ecc_lenpref.(thisfield){5}(:,1))
    % [r,p] =corr(ecc_lenpref.(thisfield){5}(:,2), ecc_lenpref.(thisfield){5}(:,1),'rows','complete');
    % title(['OPA: ', num2str(r), '/', num2str(p)]);
    % subplot(4,2,6)
    % scatter(ecc_lenpref.(thisfield){6}(:,2), ecc_lenpref.(thisfield){6}(:,1))
    % [r,p] =corr(ecc_lenpref.(thisfield){6}(:,2), ecc_lenpref.(thisfield){6}(:,1),'rows','complete');
    % title(['PPA: ', num2str(r), '/', num2str(p)]);
    % subplot(4,2,7)
    % scatter(ecc_lenpref.(thisfield){7}(:,2), ecc_lenpref.(thisfield){7}(:,1))
    % [r,p] =corr(ecc_lenpref.(thisfield){7}(:,2), ecc_lenpref.(thisfield){7}(:,1),'rows','complete');
    % title(['RSC: ', num2str(r), '/', num2str(p)]);
    % totalTitle = ['sub', num2str(isub), ' preferred length x eccentricity'];
    % sgtitle(totalTitle);

    % saveas(gcf,['/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Length/len_hist/prefLen_E', num2str(isub), '.png']);

end

%% combine all participants - preferred length from coef
% allvoxPref={};
% for visualRegion = 1:7
%     currRoi = [];
%     thisfield = combinedRoiNames{visualRegion};
%     for isub = 1:8
%         subName = ['sub', num2str(isub)];
% 
%         currRoi = [currRoi, voxPref.(subName).(thisfield)];
% 
%     end
%     allvoxPref{visualRegion} = currRoi;
% end
% 
% figure;
% subplot(4,2,1)
% histogram(allvoxPref{1});
% title('v1');
% 
% subplot(4,2,2)
% histogram(allvoxPref{2});
% title('v2');
% 
% subplot(4,2,3)
% histogram(allvoxPref{3});
% title('v3');
% 
% subplot(4,2,4)
% histogram(allvoxPref{4});
% title('v4');
% 
% subplot(4,2,5)
% histogram(allvoxPref{5});
% title('OPA');
% 
% subplot(4,2,6)
% histogram(allvoxPref{6});
% title('PPA');
% 
% subplot(4,2,7)
% histogram(allvoxPref{7});
% title('RSC');
% 
% 
% totalTitle = 'all sub number of preferred length bin';
% sgtitle(totalTitle);
% 
% saveas(gcf,'/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Length/len_hist/prefLength_allsub.png');

%% combine all participants - preffered length x eccentricity
% eccentricity < 20 %
fullecc_lenpref_all={};
for visualRegion = 1:7
    currRoi = [];
    for isub = 1:8
        subName = ['sub', num2str(isub)];
        currRoi = [currRoi; fullecc_lenpref.(subName){visualRegion}];

    end
    fullecc_lenpref_all{visualRegion} = currRoi;
end

%scatterplot
figure;
subplot(4,2,1)
scatter(fullecc_lenpref_all{1}(:,2), fullecc_lenpref_all{1}(:,1));
[r,p] =corr(fullecc_lenpref_all{1}(:,2), fullecc_lenpref_all{1}(:,1),'rows','complete');
title(['v1: ', num2str(r), '/', num2str(p)]);

subplot(4,2,2)
scatter(fullecc_lenpref_all{2}(:,2), fullecc_lenpref_all{2}(:,1));
[r,p] =corr(fullecc_lenpref_all{2}(:,2), fullecc_lenpref_all{2}(:,1),'rows','complete');
title(['v2: ', num2str(r), '/', num2str(p)]);

subplot(4,2,3)
scatter(fullecc_lenpref_all{3}(:,2), fullecc_lenpref_all{3}(:,1));
[r,p] =corr(fullecc_lenpref_all{3}(:,2), fullecc_lenpref_all{3}(:,1),'rows','complete');
title(['v3: ', num2str(r), '/', num2str(p)]);

subplot(4,2,4)
scatter(fullecc_lenpref_all{4}(:,2), fullecc_lenpref_all{4}(:,1));
[r,p] =corr(fullecc_lenpref_all{4}(:,2), fullecc_lenpref_all{4}(:,1),'rows','complete');
title(['v4: ', num2str(r), '/', num2str(p)]);

subplot(4,2,5)
scatter(fullecc_lenpref_all{5}(:,2), fullecc_lenpref_all{5}(:,1));
[r,p] =corr(fullecc_lenpref_all{5}(:,2), fullecc_lenpref_all{5}(:,1),'rows','complete');
title(['OPA: ', num2str(r), '/', num2str(p)]);

subplot(4,2,6)
scatter(fullecc_lenpref_all{6}(:,2), fullecc_lenpref_all{6}(:,1));
[r,p] =corr(fullecc_lenpref_all{6}(:,2), fullecc_lenpref_all{6}(:,1),'rows','complete');
title(['PPA: ', num2str(r), '/', num2str(p)]);

subplot(4,2,7)
scatter(fullecc_lenpref_all{7}(:,2), fullecc_lenpref_all{7}(:,1));
[r,p] =corr(fullecc_lenpref_all{7}(:,2), fullecc_lenpref_all{7}(:,1),'rows','complete');
title(['RSC: ', num2str(r), '/', num2str(p)]);

totalTitle = 'all sub preferred length x full eccentricity';
sgtitle(totalTitle);

saveas(gcf,'/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Length/len_hist/prefLen_fullE_allsub.png');

% Line plot
for visualRegion = 1:7
    V = cell(1, 8);
    for i = 1:8
        V{i} = fullecc_lenpref_all{visualRegion}(fullecc_lenpref_all{visualRegion}(:,2) == i, :);
    end

    colors = jet(8);
    figure;
    hold on;
    for i = 1:8
        histogram(V{i}(:,1), 'FaceColor', colors(i,:), 'EdgeColor', 'none', 'DisplayName', ['lengthBin', num2str(i)]);
    end
    hold off;
    legend('Location', 'best');
    title(['ecc<=20 ',combinedRoiNames{visualRegion}]);

    saveas(gcf,['/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Length/len_hist/line_prefLen_fullE_', combinedRoiNames{visualRegion}, '.png']);
end

% eccentricity <= 4.2 %
ecc_lenpref_all={};
for visualRegion = 1:7
    currRoi = [];
    for isub = 1:8
        subName = ['sub', num2str(isub)];
        currRoi = [currRoi; ecc_lenpref.(subName){visualRegion}];

    end
    ecc_lenpref_all{visualRegion} = currRoi;
end

figure;
subplot(4,2,1)
scatter(ecc_lenpref_all{1}(:,2), ecc_lenpref_all{1}(:,1));
[r,p] =corr(ecc_lenpref_all{1}(:,2), ecc_lenpref_all{1}(:,1),'rows','complete');
title(['v1: ', num2str(r), '/', num2str(p)]);

subplot(4,2,2)
scatter(ecc_lenpref_all{2}(:,2), ecc_lenpref_all{2}(:,1));
[r,p] =corr(ecc_lenpref_all{2}(:,2), ecc_lenpref_all{2}(:,1),'rows','complete');
title(['v2: ', num2str(r), '/', num2str(p)]);

subplot(4,2,3)
scatter(ecc_lenpref_all{3}(:,2), ecc_lenpref_all{3}(:,1));
[r,p] =corr(ecc_lenpref_all{3}(:,2), ecc_lenpref_all{3}(:,1),'rows','complete');
title(['v3: ', num2str(r), '/', num2str(p)]);

subplot(4,2,4)
scatter(ecc_lenpref_all{4}(:,2), ecc_lenpref_all{4}(:,1));
[r,p] =corr(ecc_lenpref_all{4}(:,2), ecc_lenpref_all{4}(:,1),'rows','complete');
title(['v4: ', num2str(r), '/', num2str(p)]);

subplot(4,2,5)
scatter(ecc_lenpref_all{5}(:,2), ecc_lenpref_all{5}(:,1));
[r,p] =corr(ecc_lenpref_all{5}(:,2), ecc_lenpref_all{5}(:,1),'rows','complete');
title(['OPA: ', num2str(r), '/', num2str(p)]);

subplot(4,2,6)
scatter(ecc_lenpref_all{6}(:,2), ecc_lenpref_all{6}(:,1));
[r,p] =corr(ecc_lenpref_all{6}(:,2), ecc_lenpref_all{6}(:,1),'rows','complete');
title(['PPA: ', num2str(r), '/', num2str(p)]);

subplot(4,2,7)
scatter(ecc_lenpref_all{7}(:,2), ecc_lenpref_all{7}(:,1));
[r,p] =corr(ecc_lenpref_all{7}(:,2), ecc_lenpref_all{7}(:,1),'rows','complete');
title(['RSC: ', num2str(r), '/', num2str(p)]);

totalTitle = 'all sub preferred length x eccentricity';
sgtitle(totalTitle);

saveas(gcf,'/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Length/len_hist/prefLen_E_allsub.png');

% Line plot
for visualRegion = 1:7
    V = cell(1, 8);
    for i = 1:8
        V{i} = ecc_lenpref_all{visualRegion}(ecc_lenpref_all{visualRegion}(:,2) == i, :);
    end

    colors = jet(8);
    figure;
    hold on;
    for i = 1:8
        histogram(V{i}(:,1), 'FaceColor', colors(i,:), 'EdgeColor', 'none', 'DisplayName', ['lengthBin', num2str(i)]);
    end
    hold off;
    legend('Location', 'best');
    title(['ecc<=4.2 ',combinedRoiNames{visualRegion}]);

    saveas(gcf,['/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Length/len_hist/line_prefLen_E_', combinedRoiNames{visualRegion}, '.png']);
end
%% 2. Grating voxel preference %%
voxModelPref = struct;
prffolder = '/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Length/prfsample_Len/';

for isub = 1:8
    subName = ['sub', num2str(isub)];
    voxModelPref.(subName)=load([prffolder 'voxModelPref_sub' num2str(isub) '.mat']);

%% preferred length from grating
    % figure;
    % subplot(4,2,1)
    % histogram(roiLen{1}(3,:));
    % title('v1');
    % 
    % subplot(4,2,2)
    % histogram(roiLen{2}(3,:));
    % title('v2');
    % 
    % subplot(4,2,3)
    % histogram(roiLen{3}(3,:));
    % title('v3');
    % 
    % subplot(4,2,4)
    % histogram(roiLen{4}(3,:));
    % title('v4');
    % 
    % subplot(4,2,5)
    % histogram(roiLen{5}(3,:));
    % title('OPA');
    % 
    % subplot(4,2,6)
    % histogram(roiLen{6}(3,:));
    % title('PPA');
    % 
    % subplot(4,2,7)
    % histogram(roiLen{7}(3,:));
    % title('RSC');
    % 
    % 
    % totalTitle = ['Grating: sub', num2str(isub), ' number of preferred length bin'];
    % sgtitle(totalTitle);
    % 
    % saveas(gcf,['/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Length/len_hist/grating_prefLength_sub', num2str(isub), '.png']);
end

%% combine all participants - preferred length from grating
allvoxModelPref={};
for roi = 1:7
    currRoi = [];
    for isub = 1:8
        subName = ['sub', num2str(isub)];

        currRoi = [currRoi, voxModelPref.(subName).roiLen{roi}(3,:)];

    end
    allvoxModelPref{roi} = currRoi;
end

figure;
subplot(4,2,1)
histogram(allvoxModelPref{1});
title('v1');

subplot(4,2,2)
histogram(allvoxModelPref{2});
title('v2');

subplot(4,2,3)
histogram(allvoxModelPref{3});
title('v3');

subplot(4,2,4)
histogram(allvoxModelPref{4});
title('v4');

subplot(4,2,5)
histogram(allvoxModelPref{5});
title('OPA');

subplot(4,2,6)
histogram(allvoxModelPref{6});
title('PPA');

subplot(4,2,7)
histogram(allvoxModelPref{7});
title('RSC');


totalTitle = 'Grating: all sub number of preferred length bin';
sgtitle(totalTitle);

saveas(gcf,'/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Length/len_hist/grating_prefLength_allsub.png');



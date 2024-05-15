% makeFigures_length.m
clear all;
cd('/home/hanseohe/Documents/GitHub/nsdOtopy/length');

prffolder = '/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Length/prfsample_Len/';  
roiNames = {'V1v','V1d','V2v','V2d','V3v','V3d','hV4','OPA','PPA','RSC'};
combinedRoiNames = {'V1','V2','V3','hV4','OPA','PPA','RSC'};

%% 1.  Coef voxel preference %%

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

%% Eccentricity scatterplot
fullecc_lenpref = struct;
ecc_lenpref = struct;
prfecc_lenpref = struct;

for isub = 1:8
    
    clearvars eccBrainbyROI prfeccBrainbyROI
    
    saveFolder = '/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Length/brainVolume/';
    load([saveFolder, 'lengthBrain_sub', num2str(isub), '.mat']);
    load([saveFolder, 'lengthBrainbyROI_sub', num2str(isub), '.mat']);

    betasfolder = ['/bwdata/NSDData/nsddata/ppdata/subj0' num2str(isub) '/func1pt8mm/'];
    eccFile = fullfile(betasfolder,'prf_eccentricity.nii.gz');
    % eccFile = fullfile(betasfolder,'roi/prf-eccrois.nii.gz');
    % 0 Unknown
    % 1 ecc0pt5
    % 2 ecc1
    % 3 ecc2
    % 4 ecc4
    % 5 ecc4+
    r2file = fullfile(betasfolder,'prf_R2.nii.gz');
    eccData = niftiread(eccFile);
    r2Data = niftiread(r2file);
    eccData(r2Data<=0) = NaN;

    eccentricity < 20
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
    % eccData(eccData > 4.2) = NaN;
    % for visualRegion = 1:7
    %     curNewBrain = newBrainbyROI(:,:,:,visualRegion);
    %     curNewBrain(curNewBrain~=-1) = eccData(curNewBrain~=-1);
    %     eccBrainbyROI(:,:,:,visualRegion) =curNewBrain;
    % end
    % for visualRegion = 1:7
    %     thisfield = ['sub', num2str(isub)];
    %     thisEcc = eccBrainbyROI(:,:,:,visualRegion);
    %     thisEcc = thisEcc(thisEcc ~= -1);
    %     thisVoxPref = newBrainbyROI(:,:,:,visualRegion);
    %     thisVoxPref = thisVoxPref(thisVoxPref ~= -1);
    %     ecc_lenpref.(thisfield){visualRegion}(:,1)=thisEcc;
    %     ecc_lenpref.(thisfield){visualRegion}(:,2)=thisVoxPref;
    % end
    
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

    % % prf-eccrois 
    % eccData(eccData == 0 | eccData == -1) = NaN;
    % for visualRegion = 1:7
    %     curNewBrain = newBrainbyROI(:,:,:,visualRegion);
    %     curNewBrain(curNewBrain~=-1) = eccData(curNewBrain~=-1);
    %     prfeccBrainbyROI(:,:,:,visualRegion) =curNewBrain;
    % end
    % for visualRegion = 1:7
    %     thisfield = ['sub', num2str(isub)];
    %     thisEcc = prfeccBrainbyROI(:,:,:,visualRegion);
    %     thisEcc = thisEcc(thisEcc ~= -1);
    %     thisVoxPref = newBrainbyROI(:,:,:,visualRegion);
    %     thisVoxPref = thisVoxPref(thisVoxPref ~= -1);
    %     prfecc_lenpref.(thisfield){visualRegion}(:,1)=thisEcc;
    %     prfecc_lenpref.(thisfield){visualRegion}(:,2)=thisVoxPref;
    % end

end

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


% prf-eccrois
prfecc_lenpref_all={};
for visualRegion = 1:7
    currRoi = [];
    for isub = 1:8
        subName = ['sub', num2str(isub)];
        currRoi = [currRoi; prfecc_lenpref.(subName){visualRegion}];

    end
    prfecc_lenpref_all{visualRegion} = currRoi;
end

%scatterplot
figure;
subplot(4,2,1)
scatter(prfecc_lenpref_all{1}(:,2), prfecc_lenpref_all{1}(:,1),'jitter','on');
[r,p] =corr(prfecc_lenpref_all{1}(:,2), prfecc_lenpref_all{1}(:,1),'rows','complete');
title(['v1: ', num2str(r), '/', num2str(p)]);

subplot(4,2,2)
scatter(prfecc_lenpref_all{2}(:,2), prfecc_lenpref_all{2}(:,1),'jitter','on');
[r,p] =corr(prfecc_lenpref_all{2}(:,2), prfecc_lenpref_all{2}(:,1),'rows','complete');
title(['v2: ', num2str(r), '/', num2str(p)]);

subplot(4,2,3)
scatter(prfecc_lenpref_all{3}(:,2), prfecc_lenpref_all{3}(:,1),'jitter','on');
[r,p] =corr(prfecc_lenpref_all{3}(:,2), prfecc_lenpref_all{3}(:,1),'rows','complete');
title(['v3: ', num2str(r), '/', num2str(p)]);

subplot(4,2,4)
scatter(prfecc_lenpref_all{4}(:,2), prfecc_lenpref_all{4}(:,1),'jitter','on');
[r,p] =corr(prfecc_lenpref_all{4}(:,2), prfecc_lenpref_all{4}(:,1),'rows','complete');
title(['v4: ', num2str(r), '/', num2str(p)]);

subplot(4,2,5)
scatter(prfecc_lenpref_all{5}(:,2), prfecc_lenpref_all{5}(:,1),'jitter','on');
[r,p] =corr(prfecc_lenpref_all{5}(:,2), prfecc_lenpref_all{5}(:,1),'rows','complete');
title(['OPA: ', num2str(r), '/', num2str(p)]);

subplot(4,2,6)
scatter(prfecc_lenpref_all{6}(:,2), prfecc_lenpref_all{6}(:,1),'jitter','on');
[r,p] =corr(prfecc_lenpref_all{6}(:,2), prfecc_lenpref_all{6}(:,1),'rows','complete');
title(['PPA: ', num2str(r), '/', num2str(p)]);

subplot(4,2,7)
scatter(prfecc_lenpref_all{7}(:,2), prfecc_lenpref_all{7}(:,1),'jitter','on');
[r,p] =corr(prfecc_lenpref_all{7}(:,2), prfecc_lenpref_all{7}(:,1),'rows','complete');
title(['RSC: ', num2str(r), '/', num2str(p)]);

totalTitle = 'all sub preferred length x prf-eccrois';
sgtitle(totalTitle);

saveas(gcf,'/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Length/len_hist/prefLen_prfE_allsub.png');

% Line plot
for visualRegion = 1:7
    V = cell(1, 8);
    for i = 1:8
        V{i} = prfecc_lenpref_all{visualRegion}(prfecc_lenpref_all{visualRegion}(:,2) == i, :);
        VV(i,:) = [sum(V{i}(:,1) == 1),sum(V{i}(:,1) == 2),sum(V{i}(:,1) == 3),sum(V{i}(:,1) == 4),sum(V{i}(:,1) == 5)];
    end

    colors = jet(8);
    figure;
    plot(VV');
    legend('lengthBin1', 'lengthBin2','lengthBin3','lengthBin4','lengthBin5','lengthBin6','lengthBin7','lengthBin8')
    legend('Location', 'best');
    title(['prf-eccrois ',combinedRoiNames{visualRegion}]);

    % saveas(gcf,['/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Length/len_hist/line_prefLen_fullE_', combinedRoiNames{visualRegion}, '.png']);
end


end


%% Size scatterplot
size_lenpref = struct;
fullsize_lenpref = struct;
prffolder = '/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Length/prfsample_Len/';
roiNames = {'V1v','V1d','V2v','V2d','V3v','V3d','hV4','OPA','PPA','RSC'};
combinedRoiNames = {'V1','V2','V3','hV4','OPA','PPA','RSC'};

for isub = 1:8
    clearvars -except isub roiNames combinedRoiNames prffolder size_lenpref fullsize_lenpref

    saveFolder = '/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Length/brainVolume/';
    load([saveFolder, 'lengthBrain_sub', num2str(isub), '.mat']);
    load([saveFolder, 'lengthBrainbyROI_sub', num2str(isub), '.mat']);

    betasfolder = ['/bwdata/NSDData/nsddata/ppdata/subj0' num2str(isub) '/func1pt8mm/'];
    sizeFile = fullfile(betasfolder,'prf_size.nii.gz');
    r2file = fullfile(betasfolder,'prf_R2.nii.gz');
    sizeData = niftiread(sizeFile);
    r2Data = niftiread(r2file);

    % size < 50 
    sizeData(sizeData > 50) = NaN;
    for visualRegion = 1:7
        curNewBrain = newBrainbyROI(:,:,:,visualRegion);
        curNewBrain(curNewBrain~=-1) = sizeData(curNewBrain~=-1);
        sizeBrainbyROI(:,:,:,visualRegion) =curNewBrain;
    end
    for visualRegion = 1:7
        thisfield = ['sub', num2str(isub)];
        thisSize = sizeBrainbyROI(:,:,:,visualRegion);
        thisSize = thisSize(thisSize ~= -1);
        thisVoxPref = newBrainbyROI(:,:,:,visualRegion);
        thisVoxPref = thisVoxPref(thisVoxPref ~= -1);
        fullsize_lenpref.(thisfield){visualRegion}(:,1)=thisSize;
        fullsize_lenpref.(thisfield){visualRegion}(:,2)=thisVoxPref;
    end


    % size <= 8.4
    clearvars sizeBrainbyROI
    sizeData(sizeData > 8.4) = NaN;
    for visualRegion = 1:7
        curNewBrain = newBrainbyROI(:,:,:,visualRegion);
        curNewBrain(curNewBrain~=-1) = sizeData(curNewBrain~=-1);
        sizeBrainbyROI(:,:,:,visualRegion) =curNewBrain;
    end
    for visualRegion = 1:7
        thisfield = ['sub', num2str(isub)];
        thisSize = sizeBrainbyROI(:,:,:,visualRegion);
        thisSize = thisSize(thisSize ~= -1);
        thisVoxPref = newBrainbyROI(:,:,:,visualRegion);
        thisVoxPref = thisVoxPref(thisVoxPref ~= -1);
        size_lenpref.(thisfield){visualRegion}(:,1)=thisSize;
        size_lenpref.(thisfield){visualRegion}(:,2)=thisVoxPref;
    end


end

%% combine all participants - preffered length x size
% size < 50
fullsize_lenpref_all={};
for visualRegion = 1:7
    currRoi = [];
    for isub = 1:8
        subName = ['sub', num2str(isub)];
        currRoi = [currRoi; fullsize_lenpref.(subName){visualRegion}];

    end
    fullsize_lenpref_all{visualRegion} = currRoi;
end

%scatterplot
figure;
subplot(4,2,1)
scatter(fullsize_lenpref_all{1}(:,2), fullsize_lenpref_all{1}(:,1));
[r,p] =corr(fullsize_lenpref_all{1}(:,2), fullsize_lenpref_all{1}(:,1),'rows','complete');
title(['v1: ', num2str(r), '/', num2str(p)]);

subplot(4,2,2)
scatter(fullsize_lenpref_all{2}(:,2), fullsize_lenpref_all{2}(:,1));
[r,p] =corr(fullsize_lenpref_all{2}(:,2), fullsize_lenpref_all{2}(:,1),'rows','complete');
title(['v2: ', num2str(r), '/', num2str(p)]);

subplot(4,2,3)
scatter(fullsize_lenpref_all{3}(:,2), fullsize_lenpref_all{3}(:,1));
[r,p] =corr(fullsize_lenpref_all{3}(:,2), fullsize_lenpref_all{3}(:,1),'rows','complete');
title(['v3: ', num2str(r), '/', num2str(p)]);

subplot(4,2,4)
scatter(fullsize_lenpref_all{4}(:,2), fullsize_lenpref_all{4}(:,1));
[r,p] =corr(fullsize_lenpref_all{4}(:,2), fullsize_lenpref_all{4}(:,1),'rows','complete');
title(['v4: ', num2str(r), '/', num2str(p)]);

subplot(4,2,5)
scatter(fullsize_lenpref_all{5}(:,2), fullsize_lenpref_all{5}(:,1));
[r,p] =corr(fullsize_lenpref_all{5}(:,2), fullsize_lenpref_all{5}(:,1),'rows','complete');
title(['OPA: ', num2str(r), '/', num2str(p)]);

subplot(4,2,6)
scatter(fullsize_lenpref_all{6}(:,2), fullsize_lenpref_all{6}(:,1));
[r,p] =corr(fullsize_lenpref_all{6}(:,2), fullsize_lenpref_all{6}(:,1),'rows','complete');
title(['PPA: ', num2str(r), '/', num2str(p)]);

subplot(4,2,7)
scatter(fullsize_lenpref_all{7}(:,2), fullsize_lenpref_all{7}(:,1));
[r,p] =corr(fullsize_lenpref_all{7}(:,2), fullsize_lenpref_all{7}(:,1),'rows','complete');
title(['RSC: ', num2str(r), '/', num2str(p)]);

totalTitle = 'all sub preferred length x size < 50';
sgtitle(totalTitle);

saveas(gcf,'/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Length/len_hist/prefLen_size50_allsub.png');

% Line plot
for visualRegion = 1:7
    V = cell(1, 8);
    for i = 1:8
        V{i} = fullsize_lenpref_all{visualRegion}(fullsize_lenpref_all{visualRegion}(:,2) == i, :);
    end

    colors = jet(8);
    figure;
    hold on;
    for i = 1:8
        histogram(V{i}(:,1), 'FaceColor', colors(i,:), 'EdgeColor', 'none', 'DisplayName', ['lengthBin', num2str(i)]);
    end
    hold off;
    legend('Location', 'best');
    title(['full size ',combinedRoiNames{visualRegion}]);

     saveas(gcf,['/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Length/len_hist/line_prefLen_size50_', combinedRoiNames{visualRegion}, '.png']);
end


% size <= 8.4 
size_lenpref_all={};
for visualRegion = 1:7
    currRoi = [];
    for isub = 1:8
        subName = ['sub', num2str(isub)];
        currRoi = [currRoi; size_lenpref.(subName){visualRegion}];

    end
    size_lenpref_all{visualRegion} = currRoi;
end

%scatterplot
figure;
subplot(4,2,1)
scatter(size_lenpref_all{1}(:,2), size_lenpref_all{1}(:,1));
[r,p] =corr(size_lenpref_all{1}(:,2), size_lenpref_all{1}(:,1),'rows','complete');
title(['v1: ', num2str(r), '/', num2str(p)]);

subplot(4,2,2)
scatter(size_lenpref_all{2}(:,2), size_lenpref_all{2}(:,1));
[r,p] =corr(size_lenpref_all{2}(:,2), size_lenpref_all{2}(:,1),'rows','complete');
title(['v2: ', num2str(r), '/', num2str(p)]);

subplot(4,2,3)
scatter(size_lenpref_all{3}(:,2), size_lenpref_all{3}(:,1));
[r,p] =corr(size_lenpref_all{3}(:,2), size_lenpref_all{3}(:,1),'rows','complete');
title(['v3: ', num2str(r), '/', num2str(p)]);

subplot(4,2,4)
scatter(size_lenpref_all{4}(:,2), size_lenpref_all{4}(:,1));
[r,p] =corr(size_lenpref_all{4}(:,2), size_lenpref_all{4}(:,1),'rows','complete');
title(['v4: ', num2str(r), '/', num2str(p)]);

subplot(4,2,5)
scatter(size_lenpref_all{5}(:,2), size_lenpref_all{5}(:,1));
[r,p] =corr(size_lenpref_all{5}(:,2), size_lenpref_all{5}(:,1),'rows','complete');
title(['OPA: ', num2str(r), '/', num2str(p)]);

subplot(4,2,6)
scatter(size_lenpref_all{6}(:,2), size_lenpref_all{6}(:,1));
[r,p] =corr(size_lenpref_all{6}(:,2), size_lenpref_all{6}(:,1),'rows','complete');
title(['PPA: ', num2str(r), '/', num2str(p)]);

subplot(4,2,7)
scatter(size_lenpref_all{7}(:,2), size_lenpref_all{7}(:,1));
[r,p] =corr(size_lenpref_all{7}(:,2), size_lenpref_all{7}(:,1),'rows','complete');
title(['RSC: ', num2str(r), '/', num2str(p)]);

totalTitle = 'all sub preferred length x size < 8.4';
sgtitle(totalTitle);

saveas(gcf,'/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Length/len_hist/prefLen_size8.4_allsub.png');

% Line plot
for visualRegion = 1:7
    V = cell(1, 8);
    for i = 1:8
        V{i} = size_lenpref_all{visualRegion}(size_lenpref_all{visualRegion}(:,2) == i, :);
    end

    colors = jet(8);
    figure;
    hold on;
    for i = 1:8
        histogram(V{i}(:,1), 'FaceColor', colors(i,:), 'EdgeColor', 'none', 'DisplayName', ['lengthBin', num2str(i)]);
    end
    hold off;
    legend('Location', 'best');
    title(['size<=8.4 ',combinedRoiNames{visualRegion}]);

    saveas(gcf,['/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Length/len_hist/line_prefLen_size8.4_', combinedRoiNames{visualRegion}, '.png']);
end



%% Eccentricity center vs periphery scatterplot
ecc_lenpref = struct;

for isub = 1:8
    
    clearvars cVSpBrainbyROI newBrainbyROIbycVSp
    
    saveFolder = '/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Length/brainVolume/';
    load([saveFolder, 'eccBrainbyROI_cVSp_sub', num2str(isub), '.mat']);
    load([saveFolder, 'lengthBrainbyROIbycVSp_sub', num2str(isub), '.mat']);


    for visualRegion = 1:14
        thisfield = ['sub', num2str(isub)];
        thisEcc = cVSpBrainbyROI(:,:,:,visualRegion);
        thisVoxPref = newBrainbyROIbycVSp(:,:,:,visualRegion);
        thisEcc = thisEcc(~isnan(thisEcc));
        thisVoxPref = thisVoxPref(~isnan(thisVoxPref));
        ecc_lenpref.(thisfield){visualRegion}(:,1)=thisEcc;
        ecc_lenpref.(thisfield){visualRegion}(:,2)=thisVoxPref;
    end

end

% combine all participants - preffered length x eccentricity

ecc_lenpref_all={};
for visualRegion = 1:14
    currRoi = [];
    for isub = 1:8
        subName = ['sub', num2str(isub)];
        currRoi = [currRoi; ecc_lenpref.(subName){visualRegion}];

    end
    ecc_lenpref_all{visualRegion} = currRoi;
end

figure;
for visualRegion = 1:14
    subplot(7, 2, visualRegion);
    scatter(ecc_lenpref_all{visualRegion}(:,2), ecc_lenpref_all{visualRegion}(:,1));
    [r, p] = corr(ecc_lenpref_all{visualRegion}(:,2), ecc_lenpref_all{visualRegion}(:,1), 'rows', 'complete');
    if mod(visualRegion, 2) == 1
        eccType = 'center';
    else
        eccType = 'periphery';
    end
    realRegion = ceil(visualRegion/2);    
    title(sprintf('%s %s: %f/%f', eccType, combinedRoiNames{realRegion}, r, p));
end
totalTitle = 'all sub preferred length x ecc by center vs. periphery';
sgtitle(totalTitle);

saveas(gcf,'/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Length/len_hist/prefLen_cVSp_allsub.png');

% Line plot
for visualRegion = 1:14
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
    if mod(visualRegion,2) == 1
        eccCond = 'center ';
    else
        eccCond = 'periphery ';
    end
    realRegion = ceil(visualRegion/2);

    title([eccCond,combinedRoiNames{realRegion}]);

    saveas(gcf,['/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Length/len_hist/line_prefLen_',eccCond,'_', combinedRoiNames{realRegion}, '.png']);


end

%% Voxel Preference center vs periphery
lenpref = struct;

for isub = 1:8
    
    clearvars  newBrainbyROIbycVSp
    
    saveFolder = '/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Length/brainVolume/';
    load([saveFolder, 'lengthBrainbyROIbycVSp_sub', num2str(isub), '.mat']);


    for visualRegion = 1:14
        subName = ['sub', num2str(isub)];
        thisVoxPref = newBrainbyROIbycVSp(:,:,:,visualRegion);
        thisVoxPref = thisVoxPref(~isnan(thisVoxPref));
        lenpref.(subName){visualRegion}=thisVoxPref;
    end

end
% combine all participants - preferred length from coef
alllenPref={};
for visualRegion = 1:14
    currRoi = [];
    for isub = 1:8
        subName = ['sub', num2str(isub)];

        currRoi = [currRoi; lenpref.(subName){visualRegion}];

    end
    alllenPref{visualRegion} = currRoi;
end

figure;
for visualRegion = 1:14
    subplot(7, 2, visualRegion);
    histogram(alllenPref{visualRegion});
    if mod(visualRegion, 2) == 1
        eccType = 'center';
    else
        eccType = 'periphery';
    end
    realRegion = ceil(visualRegion/2);
    title(sprintf('%s %s', eccType, combinedRoiNames{realRegion}));
    
end
totalTitle = 'all sub number of preferred length bin';
sgtitle(totalTitle);

saveas(gcf,'/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Length/len_hist/prefLengthbycVSp_allsub.png');


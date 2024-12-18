% getGroupMean.m
% For the group map, we computed the circular mean across subjects for each vertex in V1–V4, 
% weighted by the full model R2 values.
clear all;

filedir = '/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Orientation/surfaceData/';
condition = {'old', 'ori', 'control'};

addpath(genpath('/usr/local/freesurfer/7.4.1/matlab'));
[~,M,mr] = load_mgh('/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Orientation/surfaceData/sub1/oldBrain_sub1_lh_fsaverage.mgh');
%save_mgh(vol, 'orig.stripped.mgz', M,mr);
%result = MRIread('/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Orientation/surfaceData/sub1/angleBrain_sub1_lh_fsaverage.mgh');

% 50 percet cutoff for old R2 (participant averaged)
cutoff_lh = 0.02;
cutoff_rh = 0.02;

for curcond = 1:3
    %% get all subjects volume in radians
    % what to do with -1 and 0?
    for isub=1:8
        vol_lh(:,isub) = load_mgh([filedir, 'sub', num2str(isub), '/', condition{curcond}, 'Brain_sub', num2str(isub), '_lh_fsaverage.mgh']);
        vol_rh(:,isub) = load_mgh([filedir, 'sub', num2str(isub), '/', condition{curcond}, 'Brain_sub', num2str(isub), '_rh_fsaverage.mgh']);
    end
    
    % multiply by 2
    vol_lh = vol_lh*2;
    vol_rh = vol_rh*2;
    vol_lh_r = deg2rad(vol_lh);
    vol_rh_r = deg2rad(vol_rh);

    %% get all subjects full model R2 values
    % what to do with -1 and 0?
    for isub=1:8
        R2_lh(:,isub) = load_mgh([filedir, 'sub', num2str(isub), '/', condition{curcond}, 'BrainR2_sub', num2str(isub), '_lh_fsaverage.mgh']);
        R2_rh(:,isub) = load_mgh([filedir, 'sub', num2str(isub), '/', condition{curcond}, 'BrainR2_sub', num2str(isub), '_rh_fsaverage.mgh']);
    end


    %% weighted circ mean
    prefAngle_lh = [];
    prefAngle_rh = [];
    for ivox=1:size(vol_lh,1)
        prefAngle_lh(ivox) = circ_mean(vol_lh_r(ivox,:)',R2_lh(ivox,:)', 'omitnan');
    end
    for ivox=1:size(vol_rh,1)
        prefAngle_rh(ivox) = circ_mean(vol_rh_r(ivox,:)',R2_rh(ivox,:)', 'omitnan');
    end
   
    prefAngle_lh = mod(prefAngle_lh, 2*pi);
    prefAngle_lh = prefAngle_lh./2;%range 0 to pi.
    prefAngle_lh = prefAngle_lh/pi*180;
    prefAngle_lh = prefAngle_lh';
    prefAngle_rh = mod(prefAngle_rh, 2*pi);
    prefAngle_rh = prefAngle_rh./2;%range 0 to pi.
    prefAngle_rh = prefAngle_rh/pi*180;
    prefAngle_rh = prefAngle_rh';

    %% only top 50% full model R2 of old
    meanR2_lh = mean(R2_lh,2, "omitnan");
    meanR2_rh = mean(R2_rh,2, "omitnan");
    % cutoff_lh = prctile(meanR2_lh, 50);
    % cutoff_rh = prctile(meanR2_rh, 50);
    top50_prefAngle_lh = prefAngle_lh;
    top50_prefAngle_lh(meanR2_lh < cutoff_lh) = 0;
    top50_prefAngle_rh = prefAngle_rh;
    top50_prefAngle_rh(meanR2_rh < cutoff_rh) = 0;
    %% save
    fileName_lh = [filedir, condition{curcond}, 'Brain_groupmean_lh_fsaverage.mgh'];
    fileName_rh = [filedir, condition{curcond}, 'Brain_groupmean_rh_fsaverage.mgh'];
    save_mgh(prefAngle_lh, fileName_lh, M,mr);
    save_mgh(prefAngle_rh, fileName_rh, M,mr);

    fileName_lh = [filedir, condition{curcond}, 'Brain_groupmean_top50_lh_fsaverage.mgh'];
    fileName_rh = [filedir, condition{curcond}, 'Brain_groupmean_top50_rh_fsaverage.mgh'];
    save_mgh(top50_prefAngle_lh, fileName_lh, M,mr);
    save_mgh(top50_prefAngle_rh, fileName_rh, M,mr);
end





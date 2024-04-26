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



cd '/home/hanseohe/Documents/GitHub/nsdOtopy'

load('/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/prfsample_Ori/voxModelPref_sub1.mat');

ourBrain = visRoiData;
ourBrain(ourBrain == 2) = 1;
ourBrain(ourBrain == 3 | ourBrain == 4) = 2;
ourBrain(ourBrain == 5 | ourBrain == 6) = 3;
ourBrain(ourBrain == 7) = 4;
angleBrain = ourBrain;
angleBrain(angleBrain > 0) = 0;
angleBrain(ourBrain == 1) = roiOri{1,1}(3,:);
angleBrain(ourBrain == 2) = roiOri{1,2}(3,:);
angleBrain(ourBrain == 3) = roiOri{1,3}(3,:);
angleBrain(ourBrain == 4) = roiOri{1,4}(3,:);
angleBrain = angleBrain / pi *180;
angleBrain(angleBrain > 0) = 180 - angleBrain(angleBrain > 0);
angleBrain(angleBrain < 0) = -1;

%hist(angleBrain);


save('angleBrain_new.mat', 'angleBrain');

%%
clear all;
cd '/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/brainVolume'
load("angleBrain_new.mat");
size(angleBrain)

%3dinfo betas_session01.nii.gz 
%3dcalc -a betas_session01.nii.gz[1] -expr a -prefix oneBeta

[err,V,Info] = BrikLoad('oneBeta+orig');
Info.RootName = 'angleBrain_new+orig';
opt.Prefix = 'angleBrain_new';
WriteBrik(angleBrain,Info,opt);

%%
load('/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/prfsample_Ori/voxModelPref_sub1.mat');
roiOri_New = roiOri;
load('/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/prfsample/voxModelPref_sub1.mat');

original = roiOri{1,1}(3,:)';
new = roiOri_New{1,1}(3,:)';
[r,p] = corr(original, new);
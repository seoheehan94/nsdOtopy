%%cd '/home/hanseohe/Documents/GitHub/nsdOtopy'
cd '/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/brainVolume';
corrList = [];
for sub = 1:8


%%
fileName = ['/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/prfsample_Ori/voxModelPref_sub', num2str(sub), '.mat'];
load(fileName);

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

saveName = ['angleBrain_sub', num2str(sub), 'mat'];
save(saveName, 'angleBrain');

%%
%clear all;
%cd '/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/brainVolume'
%load("angleBrain_new.mat");
size(angleBrain)

%3dinfo betas_session01.nii.gz 
%3dcalc -a betas_session01.nii.gz[1] -expr a -prefix oneBeta
command = ['3dinfo betas_session01_sub', num2str(sub), '.nii.gz'];
system(command);
command = ['3dcalc -a betas_session01_sub', num2str(sub), '.nii.gz[1] -expr a -prefix oneBeta_sub', num2str(sub)];
system(command);

currOneBeta = ['oneBeta_sub', num2str(sub), '+orig'];
[err,V,Info] = BrikLoad(currOneBeta);
Info.RootName = ['angleBrain_sub', num2str(sub), '+orig'];
opt.Prefix = ['angleBrain_sub', num2str(sub)];
WriteBrik(angleBrain,Info,opt);

% compare original and new
roiOri_New = roiOri;
origfileName = ['/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/prfsample/voxModelPref_sub', num2str(sub), '.mat'];
load(origfileName);

original = roiOri{1,1}(3,:)';
new = roiOri_New{1,1}(3,:)';
[r,p] = corr(original, new);

corrList = [corrList; r, p];
end
%%
% load('/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/prfsample_Ori/voxModelPref_sub1.mat');
% roiOri_New = roiOri;
% load('/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/prfsample/voxModelPref_sub1.mat');
% 
% original = roiOri{1,1}(3,:)';
% new = roiOri_New{1,1}(3,:)';
% [r,p] = corr(original, new);

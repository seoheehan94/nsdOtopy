cd '/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/stimuli';
nsdfolder = '/bwdata/NSDData/nsddata/experiments/nsd/';
nsdDesignFilename = fullfile(nsdfolder, 'nsd_expdesign.mat');
nsdDesign = load(nsdDesignFilename);
% allImgs = nsdDesign.sharedix; %indices of the shared 1000 images

% mkdir subj06
% mkdir subj07
% mkdir subj08
for isub= 6:8
    fileList = '';
    allImgs = nsdDesign.subjectim(isub,nsdDesign.masterordering);%indices of all 10000 images used for this subject
    allImgs = unique(allImgs);

    for curImg = 1: length(allImgs)
        curImgName = ['pyramid/pyrImg', num2str(allImgs(curImg)), '.mat'];
        folderName = ['pyramid/pyrImg/subj0', num2str(isub)];
        movefile(curImgName, folderName)


    end
    % command = ['scp -r ', fileList, 'hanseohe@beluga4.computecanada.ca:/home/hanseohe/scratch/stimuli/pyramid'];
    % system(command);

end
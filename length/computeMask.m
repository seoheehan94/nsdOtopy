% Define image size
imageSize = 425;

% Define the desired number of pixels for each region
numPixelsPerRegion = sum(sum(true(imageSize)));

% Create meshgrid
[X, Y] = meshgrid(1:imageSize, 1:imageSize);

% Define circle parameters
centerX = imageSize / 2;
centerY = imageSize / 2;

% Initialize radius
degree = 2;
radius = 425/8.4*degree;
% radius = 169;
% radius = 170;

% Initialize masks
centerMask = false(imageSize);
peripheryMask = false(imageSize);

% Calculate the desired area for each region
desiredAreaPerRegion = numPixelsPerRegion / 2;

% Update the circular mask for the center
centerMask = (X - centerX).^2 + (Y - centerY).^2 <= radius^2;

% Calculate the area of the center region
centerArea = sum(centerMask(:));

% Calculate the area of the periphery region
peripheryArea = imageSize^2 - centerArea;

peripheryArea - centerArea

% Create periphery mask by excluding the center region
peripheryMask = ~centerMask;

% Visualize the masks
figure;
subplot(1, 2, 1);
imshow(centerMask);
title('Center Mask');

subplot(1, 2, 2);
imshow(peripheryMask);
title('Periphery Mask');

save('/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Length/Centervs.Periphery/centerMask.mat', 'centerMask');
save('/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Length/Centervs.Periphery/peripheryMask.mat', 'peripheryMask');

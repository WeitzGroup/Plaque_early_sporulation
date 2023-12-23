%This code reads image nov_3_70.png and analyzes the plaque areas in
%pixels. It saves all the plaque areas in 'NS_nov3_areas.csv'
clear
clc
close all

%% read image and convertto bw
A = imread("DSM/nov_3_70.png");
A = rgb2gray(A);
imshow(A)
%%
A_adj = imgaussfilt(A,1); % reduce noise
%histogram(A_adj)
%%
A_mask = A_adj > 160; % use general threshold to binarize 
A_mask = imfill(A_mask, 'holes');
imshowpair(A_adj, A_mask)
%% get connected components and properties
A_bw = bwlabel(A_mask); %matrix of double
A_CC = regionprops(A_mask, 'all');

%% Try watershed algorithm, oversegmentation problem
% code adapted from mathworks
D = bwdist(~A_bw);
D = -D;
D = imgaussfilt(D,3); % reduce noise

L = watershed(D);
L(~A_bw) = 0;
rgb = label2rgb(L,'jet',[.5 .5 .5]);
figure()
imshow(A)
figure()
imshow(rgb)
title('Watershed Transform')

A_mask = L > 0;
A_bw = bwlabel(A_mask); %matrix of double
CC_M_new = regionprops(A_mask, 'all');
%% select CC based on area and maybe circularity
CC_circ = [CC_M_new.Circularity];
CC_areas = [CC_M_new.Area];
% circularity larger than a threshold 
allowableCircularityIndexes = CC_circ > 0.45; 
%are between two numbers 
allowableAreaIndexes = (CC_areas >=400) & (CC_areas < 10^4); 
keeperIndexes = find(allowableAreaIndexes & allowableCircularityIndexes);
A_mask_new = ismember(A_bw, keeperIndexes); %logical matrix

% Re-label with only the CC we selected
A_bw_new = bwlabel(A_mask_new, 8);     % Label each CC 
% Plot before and after selection
A_CC_new_color = label2rgb (A_bw_new, 'hsv', 'k', 'shuffle'); % pseudo random color labels

%figure()
%imshowpair(A, A_mask)
figure()
imshowpair(A, A_mask_new)
%% Plot connected components with numbers
A_bw = bwlabel(A_mask_new); %matrix of double
CC_M_new = regionprops(A_mask_new, 'all');
n_CC_new = size(CC_M_new, 1);
figure()
imshow(A_CC_new_color)
hold on 
textFontSize = 14;	% Used to control size of "blob number" labels put atop the image.
labelShiftX = -7;	% Used to align the labels in the centers of the coins.
blob_areas = zeros(1, n_CC_new);

for k = 1 : n_CC_new           % Loop through all blobs.
	thisBlobsPixels = CC_M_new(k).PixelIdxList;  % Get list of pixels in current blob.
	meanGL = mean(A(thisBlobsPixels)); % Find mean intensity (in original image!)
	
	blobArea = CC_M_new(k).Area;		% Get area.
	blobPerimeter = CC_M_new(k).Perimeter;		% Get perimeter.
	blobCentroid = CC_M_new(k).Centroid;		% Get centroid one at a time
	blob_areas(k) = blobArea;					% Compute ECD - Equivalent Circular Diameter.
	% Put the "blob number" labels on the "boundaries" grayscale image.
	text(blobCentroid(1) + labelShiftX, blobCentroid(2), num2str(k), 'FontSize', textFontSize, 'FontWeight', 'Bold', 'Color', [1 1 1]);
end
%% Delete some plaques by hand :(
%delete 38 98 96 and 94 in decreasing order
A_bw_new = bwlabel(A_mask_new, 8);
allowableIndexes = 1: size(CC_M_new, 1);


A_mask_fin = ismember(A_bw_new, allowableIndexes);
A_bw_fin = bwlabel(A_mask_fin, 8);  
% Plot before and after selection with pseudo random color labels
A_CC_color_fin = label2rgb (A_bw_fin, 'hsv', 'k', 'shuffle'); 
A_CC_fin = regionprops(A_mask_fin, 'all');

figure()
imshowpair(A_CC_new_color, A_CC_color_fin, 'montage')
%%
figure('Position', [100 200 1200 400])
tiledlayout(1,3, 'Padding', 'none', 'TileSpacing', 'compact'); 
nexttile
imshow(A)
title('Original', 'FontSize', 16)
nexttile
imshowpair(A_adj, A_mask_fin)
title('Processed', 'FontSize', 16)
nexttile
imshow(A_adj)
hold on
centroids = cat(1,A_CC_fin.Centroid);
plot(centroids(:,1), centroids(:,2), 'rx','MarkerSize', 5, ...
    'LineWidth', 1)
title('Find and overlay centers', 'FontSize', 16)
%%
areas = cat(1,A_CC_fin.Area);
%2543 = 9cm
%%
writematrix(areas, 'NS_nov3_areas.csv')
%%
histogram(areas)
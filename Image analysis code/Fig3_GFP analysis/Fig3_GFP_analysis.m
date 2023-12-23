% This script analyzes the GPF distribution in image '30_GFP_enhanced.jpg'.
% It uses the background.png to adjust for background fluorescence and
% plots all intermediate steps as well

% Author: Andreea Magalie amagalie@gatech.edu
% Last reviewd: Dec 9th 2023
clear
clc
%% Read both images
background_GFP = imread('background.png');
A_color = imread('30_GFP_enhanced.jpg');
A = A_color(:,:,2); % take just the green component
A_len = length(background_GFP);
%% Compute background fluorescence
% We rotate the background image 4 times and  
background_all = -6*ones(4, A_len, A_len);
for x = 1:4
    temp_background = imrotate(background_GFP(:,:,2), 90*x); % rotate 90 degrees
    temp_background_mask = temp_background>55; % greate a mask to take away GFP from bacteria
    mask_uint8 = uint8(temp_background_mask);
    % delete the GFP pixels
    background_adj = temp_background - temp_background.*mask_uint8;
    % dilate the missing pixels to get a better fill later
    badPixels = imdilate(temp_background_mask, ones(15)); 
    % fill missing pixels with adjacent average 
    new_background = regionfill(background_adj, badPixels);
    %final smoothing
    background_all(x, :, :) = imgaussfilt(new_background, 2);
end
%% Average over all rotations
background_all = uint8(background_all);
fin_background = mean(background_all,1);
fin_background = uint8(reshape(fin_background, A_len, A_len));
imshow(fin_background, [10, 50])
%% Show GFP image without background
A_new = A-fin_background;
A_color(:,:,2) = A_new;
fig = figure();
brightness = 15;
imshow(A_color+brightness)
hold on
plot(150:259, ones(1,110)*1400, 'w-','LineWidth', 3)
text(130, 1450, '\boldmath{$20 \mu m$}', 'Color', 'white', 'FontSize', 18, 'Interpreter', 'latex')
set(gcf, 'Color', 'w');
%% Plot intermediate steps 
new_background = background_GFP;
new_background(:,:,2) = fin_background;
figure('Position', [0 300 800 800])
subplot(2,2,1)
imshow(background_GFP)
title('Original background')
subplot(2,2,2)
imshow(new_background)
title('Adjusted background')
subplot(2,2,3)
imshow(A_color+brightness)
title('Original image')
subplot(2,2,4)
imshow(A_color+brightness)
title('Adjusted image')
%% Blurr new GFP image
window_size = floor(A_len/25);
figure()
subplot(1,2,1);
imshow(A_new, [0 217]);
title('Original Image', 'FontSize', 15);
A_blurr = imfilter(A_new, ones(window_size)/window_size^2, 'symmetric');
subplot(1,2,2);
imshow(A_blurr, [0, 217]);
title('Blurred Image', 'FontSize', 15);
%% Find peaks of plaque edge along x and y axis
% Not the prettiest solution but it works 
x_sum = sum(A_blurr, 2); % col/ row sum of blurred image
y_sum = sum(A_blurr, 1);

interval_1 = 1:700; % general interval of peak 1
interval_2 = 700:1200; % general interval of peak 2

x_max_1_val = max(x_sum(interval_1));
x_max_2_val = max(x_sum(interval_2));
x_peak_1 = find(x_sum(interval_1) == x_max_1_val);
x_peak_2 = find(x_sum(interval_2) == x_max_2_val);

y_max_1_val = max(y_sum(interval_1));
y_max_2_val = max(y_sum(interval_2));
y_peak_1 = find(y_sum(interval_1) == y_max_1_val);
y_peak_2 = find(y_sum(interval_2) == y_max_2_val);

% find center at average of peaks
center_p = [x_peak_1 + 700+x_peak_2, y_peak_1 + 700+y_peak_2];
center_p = center_p/2;
%% Plot summary of center and x and y sums 
figure('Position', [100 300 1000 400])
y_tics = linspace(0, 5*10^4,100);
subplot(1,2,2)
plot(x_sum, 'LineWidth', 3)
hold on
plot(y_sum, 'LineWidth', 3)
plot(x_peak_1, x_max_1_val, 'o', 'LineWidth', 3, 'Color', ...
    [0 0.4470 0.7410] ,'MarkerSize', 16)
plot(x_peak_2+700, x_max_2_val, 'o', 'LineWidth', 3, 'Color', ...
    [0 0.4470 0.7410] ,'MarkerSize', 16)
plot(y_peak_1, y_max_1_val, 'o', 'LineWidth', 3, 'Color', ...
    [0.8500 0.3250 0.0980] ,'MarkerSize', 16)
plot(y_peak_2+700, y_max_2_val, 'o', 'LineWidth', 3, 'Color', ...
    [0.8500 0.3250 0.0980] ,'MarkerSize', 16)
xlabel('Location (pixels)')
ylabel('GFP intensity')
set(gca, 'FontSize', 16)
legend('x axis GFP sum', 'y axis GFP sum')
xlim([0 1608])
subplot(1,2,1)
A_blurr_all = imfilter(A_color+brightness, ones(window_size)/window_size^2, 'symmetric');
A_new = A_blurr_all;
imshow(A_new)
hold on
title('Blurred image')
plot(center_p(2), center_p(1), 'wx','LineWidth', 2.5, 'Markersize', 12)
set(gca, 'FontSize', 16)
%% compute distance from center of plaque to every other pixel
dist_tocenter = zeros(A_len, A_len);
for x = 1:A_len
    for y = 1:A_len
        dist_tocenter(x, y) = sqrt((x-center_p(1))^2 + (y-center_p(2))^2);
    end
end
%% divide pixels based on their distance to center in a number of bins 
A_blurr_copy = im2double(A_blurr);
min_dist = min(min(dist_tocenter));
max_dist = max(max(dist_tocenter));
n_bins = 400; % numbers of bins to divide distance to center
int_values = linspace(min_dist, max_dist, n_bins);

gfp_plot = NaN(n_bins,3e5);
count_plot = zeros(n_bins,1);
A_partition = zeros(A_len, A_len);

for i_bin = 1 : n_bins-1
     % find all pixels st distance to center is in a given bin
    [x_vals, y_vals] = find(dist_tocenter < int_values(i_bin+1) & ...
    dist_tocenter > int_values(i_bin));

    for mean_GFP = 1: length(x_vals)
        gfp_plot(i_bin, mean_GFP) = A_blurr_copy(x_vals(mean_GFP), y_vals(mean_GFP)); % value of GFP in bin
        count_plot(i_bin) = count_plot(i_bin) + 1; % how many points in bin
        A_partition(x_vals(mean_GFP), y_vals(mean_GFP)) = i_bin; % partition by bin number 
    end
end
A_partition = uint8(A_partition);
%% Plot
mean_GFP = nanmean(gfp_plot, 2);
err = std(gfp_plot, 0, 2, 'omitnan');

curve1 = mean_GFP + err;
curve2 = mean_GFP - err;

curve1 = curve1(1:end-1);
curve2 = curve2(1:end-1);
mean_GFP = mean_GFP(1:end-1);
int_values = linspace(min_dist, max_dist, n_bins);
int_values = int_values(1:end-1);
figure('Position', [100 300 460 460])
hold on
p1 = patch([int_values fliplr(int_values)], [curve1' fliplr(curve2')], ...
    [0.4160 0.7540 0.1880],'HandleVisibility', 'off');
p1.FaceAlpha = .2;
p1.EdgeColor = 'none';
grid on

plot(int_values, mean_GFP, 'LineWidth', 3, 'Color', [0.4160 0.7540 0.1880])
plot(int_values, curve1,'--', 'Color',[.3 .7 .1], 'LineWidth',2)
plot(int_values, curve2, '--g', 'LineWidth',2, 'Color',[.3 .7 .1]);
set(gca, 'FontSize', 18)
legend('Mean', 'Standard deviation')
xlabel('Distance from center (\mu m)')
ylabel('GFP intensity')
ylim([0 0.21])
xlim([0 1100])
%294.87microns  = 1608 pixels
xlabels = [0, 25, 75, 125, 175];
xind = xlabels * 1608/295;
xticks(xind)
xticklabels({xlabels})

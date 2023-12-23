% This script uses the plaque centers to define a search window and compute
% plaque areas at earlier time points. It accounts for camera movement
% using the edges_nov19.mat dataset
clc
clear
close all
%% Load centers of plaques at end of time lapse
stats = readmatrix('Data/WT_nov19_centers.csv');
n_plaques = length(stats); % number of plaques
start_point = 1; % start frame of plaque size analysis 
fin_point = 179; % end frame of plaque size analysis 
len_frames = length(start_point:fin_point);
areas = zeros(len_frames, n_plaques);
n_CC = zeros(len_frames, n_plaques);

load 'Data/edges_nov19.mat'
end_center = edges(:, :,fin_point); % edges at the last frame
%% For each frame before end, go through each plaque center and compute area 
for n_i = start_point:fin_point
    % create string name of file that needs to be opened
    if n_i < 10
        mystring = strcat('cropped/WT_DSM_SPO1_00',num2str(n_i),'.jpg');
    elseif n_i < 100
        mystring = strcat('cropped/WT_DSM_SPO1_0',num2str(n_i),'.jpg');
    else
        mystring = strcat('cropped/WT_DSM_SPO1_',num2str(n_i),'.jpg');
    end
    A = imread(mystring);
    se = strel('disk',25);
    Aprime = imgaussfilt(A,1); % reduce noise
    Aprime = imopen(Aprime, se);
    A_adj = imadjust(A-Aprime);
    A_adj = imgaussfilt(A_adj,1); % reduce noise    
    %imshow(A_adj)
    drift_y = edges(1,:, n_i) - end_center(1, :); %drift in y direction
    drift_x = edges(2,:, n_i) - end_center(2, :); %drift in x direction
    drift_x = floor(mean(drift_x));
    drift_y = floor(mean(drift_y));
    for i = 1 : n_plaques
        l = stats(i, 2);
        L = stats(i, 1); % x dim
        c_x = stats(i, 3) + drift_x;  %xdim
        c_y = stats(i, 4) + drift_y; %ydim
        rect = [(c_x-L/2), (c_y-l/2), L, l];
        temp_crop = imcrop(A_adj, rect);
        pl_area = analyze_single_pl_wt(temp_crop);
        areas(n_i-start_point+1, i) = sum(sum(pl_area));
        CC = bwconncomp(pl_area);
        n_CC(n_i-start_point+1, i) = CC.NumObjects;

        if n_i<=33 || CC.NumObjects==0
            areas(n_i-start_point+1, i) = 0;
        elseif CC.NumObjects>1 && n_i<70
            areas(n_i-start_point+1, i) = 0;
        else
            areas(n_i-start_point+1, i) = max(cellfun(@numel,CC.PixelIdxList));
        end
    end
end
%% Show mean and individual trajectories
figure('Position', [100 300 1200 1100])
subplot(2,2,1)
plot(start_point:fin_point, mean(areas, 2), 'LineWidth', 2)
set(gca, 'FontSize', 16)
ylabel('Mean plaque size')
xlabel('Time')

subplot(2,2,2)
plot(start_point:fin_point, areas(:,:), 'LineWidth', 1)
set(gca, 'FontSize', 16)
ylabel('Plaque size')
xlabel('Time')

subplot(2,2,3)
plot(start_point:fin_point, mean(n_CC, 2), 'LineWidth', 2)
set(gca, 'FontSize', 16)
ylabel('Mean number of CC')
xlabel('Time')

subplot(2,2,4)
plot(start_point:fin_point, n_CC(:,:), 'LineWidth', 1)
set(gca, 'FontSize', 16)
ylabel('Number of CC')
xlabel('Time')
%% save data
%save('Data/WT_areas_v2.mat', 'areas')
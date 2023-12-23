clc
clear
% this code loads the plaque data from NS_nov3_areas.csv and
% S_nov3_areas.csv and plots it using the imported function violinplot
%% Load areas in pixels from previous analysis
load 'NS_nov3_areas.csv'
NS_areas = NS_nov3_areas; %2335 9cm
load 'S_nov3_areas.csv'
S_areas = S_nov3_areas; %2543 9cm
%% transform from area in pixels to radius in mm
% From inspection of original image, we have 2335 pixels for the diameter
% of the petri dish for NS and 2543 pixels for S
NS_areas = NS_areas *(90/2335)^2; % 1 pixel = 90/2335 mm
NS_radius = sqrt(NS_areas/pi); % get radius
S_areas = S_areas *(90/2543)^2;
S_radius = sqrt(S_areas/pi);
%%
mean_NS = mean(NS_radius);
mean_S = mean(S_radius);
%% Prep datasets for violin plot
tot_areas = [S_areas; NS_areas];
tot_radius = sqrt(tot_areas/pi);
radius_labels = cell(length(tot_radius),1);
for i = 1:length(S_areas)
    radius_labels(i)={'S'};
end
for i = 1:length(NS_areas)
    radius_labels(length(S_areas)+i)={'NS'};
end
%% Use violin plot
figure('Position', [100 300 550 460])
hold on
t_S = 0.7:0.01:1.3;
t_NS = 1.72:0.01:2.28;
%1.446 and 0.6762 are the median values from the plot
plot(t_S,  1.446* ones(1, length(t_S)), 'LineWidth', 3.5, 'Color', [0, 0.4470, 0.7410])
plot(t_NS,  0.6762* ones(1, length(t_NS)), 'LineWidth', 3.5, 'Color', ...
    [0.8500, 0.3250, 0.0980])
violinplot(tot_radius,radius_labels, 'EdgeColor', [.6 .6 .6]...
    , 'ViolinAlpha', .4, 'ShowMedian', true, 'ShowBox', false, ...
    'ShowWhiskers', false)

legend('SPOIIE median', 'WT median')
ylabel('Plauqe radius (mm)')
set(gca, 'FontSize', 18)
grid on
%% Ttest to check that distributions have different means 
S_radius = sqrt(S_areas/pi);
NS_radius = sqrt(NS_areas/pi);

[h, p] = ttest2(S_radius, NS_radius);

% This script plots the final plaque areas over time
clc
clear
close all
%% load data
start_point = 1; 
fin_point = 179;
window_size = 7;
pixel_size = 90/1300; % 90 mm is equal to 1300 pixels 
load 'Data/delta_areas_v2.mat' % loads as areas 
Delta_areas = areas * pixel_size^2; % now expressed in mm^2
Delta_areas_adjusted = movmean(Delta_areas, window_size); % use a moving mean to reduce noise
SPOIIE_radius = sqrt(Delta_areas_adjusted/pi); % from area to radius
SPOIIE_mean_radius = mean(SPOIIE_radius,2); % mean radius in mm

load 'Data/WT_areas_v2.mat'
WT_areas = areas * pixel_size^2;
WT_areas_adjusted = movmean(WT_areas, window_size);
WT_radius = sqrt(WT_areas_adjusted/pi);
WT_mean_radius = mean(WT_radius, 2);

%%
% there are 179 frames, taken every 5 minutes. It is assumed that frame 1 is
% taken at minute 30 after experiment has started (time from when bacteria 
% and viruses mix until the camera takes the first shot)
t_vals = ((start_point:fin_point-1)*5 + 30)/60; % time in hours
figure('Position', [100 300 550 460])
hold on
err = std(SPOIIE_radius, 0, 2);
y1 = SPOIIE_mean_radius + err;
y2 = SPOIIE_mean_radius - err;
p1 = patch([t_vals fliplr(t_vals)],[y1' fliplr(y2')], [0, 0.4470, 0.7410],...
    'HandleVisibility', 'off');
p1.FaceAlpha = .2; 
p1.EdgeColor = 'none' ; 

plot(t_vals, SPOIIE_mean_radius + err, '--','LineWidth', 2, 'Color',...
    [0, 0.4470, 0.7410 .5], 'HandleVisibility', 'off')
plot(t_vals, SPOIIE_mean_radius - err, '--','LineWidth', 2, 'Color',...
    [0, 0.4470, 0.7410 .5], 'HandleVisibility', 'off')

err_S = std(WT_radius, 0, 2);
plot(t_vals,WT_mean_radius+ err_S, '--',...
    'LineWidth', 2, 'Color', [0.8500, 0.3250, 0.0980 .5], 'HandleVisibility', 'off')
plot(t_vals,WT_mean_radius- err_S, '--',...
    'LineWidth', 2, 'Color', [0.8500, 0.3250, 0.0980 .5], 'HandleVisibility', 'off')
set(gca, 'FontSize', 18)
ylabel('Mean plaque radius (mm)')
xlim([t_vals(1) t_vals(179)])
ylim([0 1.8])
xlabel('Time (hours)')
xticks([1 3 5 7 9 11 13 15])

grid on
legend('Location', 'NW')
y1 = WT_mean_radius + err_S;
y2 = WT_mean_radius - err_S;
p1 = patch([t_vals fliplr(t_vals)],[y1' fliplr(y2')], [0.8500, 0.3250, 0.0980],...
    'HandleVisibility', 'off');
p1.FaceAlpha = .2 ; 
p1.EdgeColor = 'none' ; 
plot(t_vals, SPOIIE_mean_radius,'Color',[0, 0.4470, 0.7410], 'LineWidth',...
    3, 'DisplayName','SPOIIE');
plot(t_vals, WT_mean_radius, 'LineWidth', 3,...
    'Color', [0.8500, 0.3250, 0.0980],'DisplayName','WT');

xticks([1 3 5 7 9 11 13 15])
yticks([0, 0.4, 0.7, 1, 1.3, 1.6 ])
%% save data in mm
%save('Data/SPOIIE_mean_radius.mat', 'SPOIIE_mean_radius')
%save('Data/WT_mean_radius', 'WT_mean_radius')
%% Plot individual trajectories
figure()
plot(t_vals, SPOIIE_radius, ...
    'LineWidth', 1, 'Color', [0, 0.4470, 0.7410, .2])
hold on 
plot(t_vals, WT_radius, ...
    'LineWidth', 1, 'Color', [0.85, 0.325, 0.098, .2])
set(gca, 'FontSize', 16)
%ylabel('Mean plaque radius (pixels)')
xlim([t_vals(1) t_vals(179)])
ylabel('Plaque radius (pixels)')
ylim([0 1.8])
xticks([1 3 5 7 9 11 13 15])
xlabel('Time (hours)')
grid on
title('Individual runs to be fitted')
%% fit each plaque growth
SPOIIE_slopes = ones(1, 46);
figure()
hold on
% for each plaque based on the end value find the 50-70% interval and fit a
% line to that
for i=1:46
ind1 = find(SPOIIE_radius(:,i)>0.5*SPOIIE_radius(end,i), 1);
ind2 = find(SPOIIE_radius(:,i)<0.7*SPOIIE_radius(end,i), 1, 'last');
ind = ind1:ind2;
temp = polyfit(t_vals(ind), SPOIIE_radius(ind,i), 1);
SPOIIE_slopes(i) = temp(1);
plot(t_vals(ind), SPOIIE_radius(ind,i))
end
%%
WT_slopes = ones(1, 100);
cmap = parula;
figure()
hold on
for i=1:100
    ind1 = find(WT_radius(:,i)>0.3*WT_radius(end,i), 1);
    ind2 = find(WT_radius(:,i)<0.7*WT_radius(end,i), 1, 'last');
    ind = ind1:ind2;
    temp = polyfit(t_vals(ind), WT_radius(ind,i), 1);
    WT_slopes(i) = temp(1);
    plot(t_vals(ind), WT_radius(ind,i), 'LineWidth', 1.5, 'Color', cmap(2*i, :))
    plot(t_vals(ind), t_vals(ind) * temp (1) + temp(2), '--', ...
        'LineWidth', 1.5,'Color', cmap(2*i, :))
end
set(gca, 'FontSize', 16)
xlabel('Time (hr)')
ylabel('Plaque size (mm)')
%% plot WT and SPOIIE slopes
histogram(WT_slopes)
hold on
histogram(SPOIIE_slopes)
set(gca, 'FontSize', 16)
xlabel('slope value')
legend('WT', 'SPOIIE')
%% boot strap slope values
n_bootstraps = 1e5;
WT_slopes_mu = WT_slopes * 10^3; % convert from mm to microm
SPOIIE_slopes_mu = SPOIIE_slopes * 10^3; % convert from mm to microm

all_slopes = [WT_slopes_mu, SPOIIE_slopes_mu]; 
n_slopes = length(all_slopes); % total number of slopes/ plaques
n_WT = length(WT_slopes_mu); % number of WT plaques/ slopes
n_SPOIIE = length(SPOIIE_slopes_mu); % number of SPOIIE plaques/ slopes
X = -1*ones(n_bootstraps,1);

for i=1:n_bootstraps
    ind_WT = randperm(n_slopes, n_WT);% select n_WT random indices
    ind_SPOIIE = setdiff(1:n_slopes, ind_WT); % get the complement indices
    A_mean = mean(all_slopes(ind_WT)); % mean of A group
    B_mean = mean(all_slopes(ind_SPOIIE)); % mean of B group
    X(i) = abs(A_mean - B_mean);
end

true_X = abs(mean(WT_slopes_mu) - mean(SPOIIE_slopes_mu));
p_val = sum(X > true_X)/n_bootstraps;
disp(p_val)

figure()
histogram(X, 'HandleVisibility', 'off')
hold on
plot(true_X * ones(1,200), linspace(0,3250, 200), '--', 'LineWidth', 3)
ylim([0 3250])
xlim([0 60])
set(gca, 'FontSize', 16)
legend('Observed difference')
xlabel('Absolute difference in mean slope ($\mu m /hr$)', 'Interpreter', 'latex')
ylabel('Counts (100,000 total)', 'Interpreter', 'latex')
%%
disp(['WT mean is ', int2str(round(mean(WT_slopes_mu)))])
disp(['WT std is ', num2str(round(std(WT_slopes_mu)))])
disp(['SPOIIE mean is ', num2str(round(mean(SPOIIE_slopes_mu)))])
disp(['SPOIIE std is ', num2str(round(std(SPOIIE_slopes_mu)))])
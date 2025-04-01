% only vertical - mrr

% Matlab script with a least-squares algorithm to carry out a vertical-component inversion.
% GEOL3048 Seismology II Practical, modified to use only the vertical sensitivity kernel.

close all
clear all

%% Load Data
% Load files with the vertical sensitivity kernel (mrr) and the real (vertical) waveform data.
folder_name = 'data';
cd(folder_name)
% Only load mrr, and the real and synthetic data if needed
filenames = {'mrr', 'real', 'synth'};
for n = 1:length(filenames)
    load(filenames{n})
end
cd ..

%% Inversion using only the vertical component
% Select input waveform data (assumed to be the vertical component)
inp_waves = real;
% inp_waves=synth;

% Build the design matrix using only the vertical sensitivity kernel.
% G is now a column vector.
G = mrr;

% Perform the least-squares inversion.
% This inverts for a single scaling parameter (a scalar) rather than a full moment tensor.
moment_scalar = (G' * G) \ (G' * inp_waves) * 1e26;

% Calculate the best-fit synthetic using only the vertical kernel.
best_fit = mrr * moment_scalar / 1e26;

%% Plotting Station Distribution
% Load station file
fileID = fopen([folder_name '/stations']);
C = textscan(fileID, '%s %s %f %f');
station = C{1};
channel = C{2};
azimuth = C{3};
distance = C{4};

% Polar plot of the source and stations
figure(1)
h = polar(deg2rad(azimuth), distance, '^');
hold on
set(h, 'markersize', 9, 'MarkerFaceColor', 'b')
xlabel('Distance from Epicenter in Degree')
ylabel('Azimuth from North in Degree')
title('Station Distribution')

view(90, -90)
h = polar(0, 0, 'pr');
set(h, 'markersize', 15, 'MarkerFaceColor', 'r')

offset = 3;
for i = 1:length(azimuth)
    angle_rad = deg2rad(azimuth(i));       % convert to radians
    [x, y] = pol2cart(angle_rad, distance(i) + offset); % convert polar to Cartesian
    % Create an outline for better visibility
    outlineOffsets = [-0.2, 0, 0.2];
    for dx = outlineOffsets
        for dy = outlineOffsets
            if ~(dx == 0 && dy == 0)
                text(x + dx, y + dy, station{i}, 'Color', 'white', 'FontSize', 8);
            end
        end
    end
    text(x, y, station{i}, 'Color', 'black', 'FontSize', 8);
end
hold off

%% Plotting the Comparison Between Input Data and Inversion Results
load([folder_name '/cutpoints'])

figure(2)
[az, az_index] = sort(azimuth); % Plot with ascending azimuth
for m = 1:length(station)
    k = az_index(m);
    subplot(7, 3, m)
    window = cutpoints(k):cutpoints(k+1);
    time = 0:(cutpoints(k+1)-cutpoints(k));
    plot(time, inp_waves(window), 'k', time, best_fit(window), 'r')
    xlim([0 800])
    ylim([-abs(max(inp_waves)) abs(max(inp_waves))]) % Uniform scale for all traces
    title([station{k} ' ' channel{k}(3) ' dist ' num2str(distance(k)) ' az ' num2str(azimuth(k))])
end

sgtitle('Waveform Comparison for All Stations')
legend('Real Data', 'Vertical Component Inversion')

%% Additional Calculations
% Since only one component is used, many advanced analyses (e.g. double-couple determination,
% eigenvalue decomposition of a full moment tensor, seismic moment calculation from six components)
% are not applicable. If needed, you could compute a simple seismic moment estimate from the scalar.
%
% For example, a simplified scalar moment (assuming mrr directly scales the vertical motion):
%
M0 = abs(moment_scalar);  
disp(['Estimated scalar seismic moment = ', num2str(M0)]);

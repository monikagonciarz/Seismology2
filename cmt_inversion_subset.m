% Matlab script with a least-squares algorithm to carry out a CMT source inversion.
% GEOL3048 Seismology II Practical week 9, Michael Frietsch and Ana Ferreira, March 2015

close all
clear all

%% Load Data Files
% Load files with six sensitivity kernels, synthetic seismograms and real data in the data/ directory
folder_name = 'data';
cd(folder_name)
filenames = {'mrr', 'mtt', 'mpp', 'mrt', 'mrp', 'mtp', 'real', 'synth'};
for n = 1:length(filenames)
    load(filenames{n})
end
cd ..

% Decide input data (synthetics or real data) 
inp_waves = real;
% inp_waves = synth;

% Build G matrix with the sensitivity kernels (sensitivity in u wrt each parameter)
G = [mrr mtt mpp mrt mrp mtp];

%% Load Station Information
% Load station file (station name, channel, azimuth and distance)
fileID = fopen([folder_name '/stations']);
C = textscan(fileID, '%s %s %f %f');
station = C{1};
channel = C{2};
azimuth = C{3};
distance = C{4};

%% Plotting Station Distribution
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
    % Create an outline by plotting the text several times in white, each slightly offset
    outlineOffsets = [-0.2, 0, 0.2];
    for dx = outlineOffsets
        for dy = outlineOffsets
            if ~(dx==0 && dy==0)
                text(x+dx, y+dy, station{i}, 'Color', 'white', 'FontSize', 8);
            end
        end
    end
    text(x, y, station{i}, 'Color', 'black', 'FontSize', 8);
end
hold off

%% Load Cutpoints
% cutpoints is assumed to have length = (number of stations + 1)
load([folder_name '/cutpoints'])

%% Filter Data for Selected Stations (for inversion and waveform comparisons)
% Define the list of selected stations (use exactly the names you want)
selectedStations = {'SJG', 'SACV', 'SFJD', 'ESK', 'BFO', 'PAB'};

% Find indices of stations that match the selected list
selIdx = find(ismember(station, selectedStations));
if isempty(selIdx)
    warning('None of the specified station names were found. Using all stations.');
    selIdx = 1:length(station);
end

% Preallocate empty arrays for the filtered data
G_selected = [];
inp_waves_selected = [];
% Also build new cutpoints for the selected stations in the filtered data
sel_cutpoints = zeros(length(selIdx)+1, 1);
sel_cutpoints(1) = 1;

for i = 1:length(selIdx)
    idx = selIdx(i);
    % Determine the window for this station from the full cutpoints
    window = cutpoints(idx):cutpoints(idx+1)-1;  % subtract 1 so that windows don’t overlap
    % Append data from G and inp_waves
    G_selected = [G_selected; G(window, :)];
    inp_waves_selected = [inp_waves_selected; inp_waves(window)];
    % Update the new cutpoints for the filtered data
    n_samples = length(window);
    sel_cutpoints(i+1) = sel_cutpoints(i) + n_samples;
end

% Ensure waveform data is a column vector
inp_waves_selected = inp_waves_selected(:);

%% Perform the Inversion with Selected Stations
% Using the standard least squares inversion: moment_tensor = (G'G)^(-1)*G'*inp_waves
moment_tensor_selected = (G_selected' * G_selected) \ (G_selected' * inp_waves_selected) * 1e26;
% Compute the best-fit synthetic waveforms for the selected inversion
best_fit_selected = G_selected * moment_tensor_selected / 1e26;

%% Find the Best Double-Couple Fault Geometry (using selected inversion)
% (Assuming mt2sdr and auxplane are available functions)
mainplane = mt2sdr(moment_tensor_selected')  % Note: transposed if needed by mt2sdr
aux_plane = auxplane(mainplane)

%% Waveform Comparison: Ordered by Distance
% For plotting, gather station information for the selected stations
sel_station = station(selIdx);
sel_channel = channel(selIdx);
sel_distance = distance(selIdx);
sel_azimuth = azimuth(selIdx);

% Create a cell array containing the indices (windows) for each selected station
win_sel = cell(length(selIdx), 1);
for i = 1:length(selIdx)
    win_sel{i} = sel_cutpoints(i):(sel_cutpoints(i+1)-1);
end

% Order by distance
[sortedD, orderD] = sort(sel_distance);
nSelected = length(orderD);
nCols = 3;
nRows = ceil(nSelected / nCols);

figure(2)
for m = 1:nSelected
    i = orderD(m);
    window = win_sel{i};
    time = 0:length(window)-1;
    % Calculate trace duration (number of samples)
    trace_duration = length(window);
    fprintf('Trace duration for station %s: %.2f samples\n', sel_station{i}, trace_duration);
    
    subplot(nRows, nCols, m)
    plot(time, inp_waves_selected(window), 'k', time, best_fit_selected(window), 'r')
    xlim([0 max(time)])
    ylim([-max(abs(inp_waves_selected)) max(abs(inp_waves_selected))])
    title([sel_station{i} ' ' sel_channel{i}(3) ' dist ' num2str(sel_distance(i)) ' az ' num2str(sel_azimuth(i))])
end
sgtitle('Waveform Comparison for Selected Stations (Ordered by Distance)')
legend('Input Data', 'Inversion')

%% Waveform Comparison: Ordered by Azimuth
[sortedAz, orderAz] = sort(sel_azimuth);
nSelected = length(orderAz);
nCols = 3;
nRows = ceil(nSelected / nCols);

figure(3)
for m = 1:nSelected
    i = orderAz(m);
    window = win_sel{i};
    time = 0:length(window)-1;
    % Calculate trace duration (number of samples)
    trace_duration = length(window);
    fprintf('Trace duration for station %s: %.2f samples\n', sel_station{i}, trace_duration);
    
    subplot(nRows, nCols, m)
    plot(time, inp_waves_selected(window), 'k', time, best_fit_selected(window), 'r')
    xlim([0 max(time)])
    ylim([-max(abs(inp_waves_selected)) max(abs(inp_waves_selected))])
    title([sel_station{i} ' ' sel_channel{i}(3) ' dist ' num2str(sel_distance(i)) ' az ' num2str(sel_azimuth(i))])
end
sgtitle('Waveform Comparison for Selected Stations (Ordered by Azimuth)')
legend('Input Data', 'Inversion')

%% Plot the Focal Mechanism Based on Selected Stations' Inversion
figure(4)
plotmt(1, 1, moment_tensor_selected')   % Transpose if plotmt expects a row vector
colormap([1 1 1; 1 0 0])  % Red for compressional, white for dilatational
title('The Focal Mechanism for Selected Stations')

%%
% Define your fault parameters for the main and auxiliary planes
mainplane_strike = 67.3381;
mainplane_dip    = 21.9487;
mainplane_rake   = -165.1510;

auxplane_strike  = 323.5226;
auxplane_dip     = 84.5033;
auxplane_rake    = -68.7178;

% Convert the fault parameters to moment tensors (1x6 vectors)
plane_mt  = sdr2mt(auxplane_strike, auxplane_dip, auxplane_rake);

% Now plot the focal mechanisms using plotmt
figure;
plotmt(1, 1, plane_mt);
colormap([1 1 1; 1 0 0]);
title('The Focal Mechanism (Station Subset)')




%% Plot an Example of Moment Tensor Component Sensitivity Kernels
figure(5)
subplot(2,3,1);
plot(mrr)
title('mrr')
xlabel('x index')
ylabel('mrr component values')
subplot(2,3,2);
plot(mpp)
title('mpp')
xlabel('x index')
ylabel('mpp component values')
subplot(2,3,3);
plot(mtt)
title('mtt')
xlabel('x index')
ylabel('mtt component values')
subplot(2,3,4);
plot(mrp)
title('mrp')
xlabel('x index')
ylabel('mrp component values')
subplot(2,3,5);
plot(mrt)
title('mrt')
xlabel('x index')
ylabel('mrt component values')
subplot(2,3,6);
plot(mtp)
title('mtp')
xlabel('x index')
ylabel('mtp component values')
sgtitle('Moment Tensor Component Sensitivity Kernels')

%% Plotting Sensitivity with Respect to a Component of the Moment Tensor
figure(6)
[~, az_index] = sort(azimuth); % Sorting all stations by ascending azimuth (using full set here)
for m = 1:length(station)
    k = az_index(m);
    subplot(7,3,m)
    window = cutpoints(k):cutpoints(k+1)-1;
    time = 0:length(window)-1;
    plot(time, mtp(window), 'b')
    xlim([0 max(time)])
    ylim([-max(abs(mtp)) max(abs(mtp))])
    title([station{k} ' ' channel{k}(3) ' dist ' num2str(distance(k)) ' az ' num2str(azimuth(k))])
end
sgtitle('Sensitivity with Respect to Component of Moment Tensor')

%% Calculating the Trace of the Moment Tensor (from the selected inversion)
% For a moment tensor vector [mrr mtt mpp mrt mrp mtp], the trace is:
trace_selected = moment_tensor_selected(1) + moment_tensor_selected(2) + moment_tensor_selected(3);
disp(['Trace of the Moment Tensor = ' num2str(trace_selected)])

%% Eigenvalue-Eigenvector Decomposition of the Seismic Moment Tensor (selected inversion)
M_sel = [moment_tensor_selected(1) moment_tensor_selected(4) moment_tensor_selected(5);
         moment_tensor_selected(4) moment_tensor_selected(2) moment_tensor_selected(6);
         moment_tensor_selected(5) moment_tensor_selected(6) moment_tensor_selected(3)];
eig(M_sel)

%% Calculating the Percentage Deviation from a Pure Double Couple Component
[V, D] = eig(M_sel);
epsilon = 100 * min(abs(diag(D))) / max(abs(diag(D)));
disp(['Percentage deviation from pure double couple = ' num2str(epsilon) ' %'])

%% Calculating the Seismic Moment (from the selected inversion)
M_rr = moment_tensor_selected(1);
M_tt = moment_tensor_selected(2);
M_pp = moment_tensor_selected(3);
M_rt = moment_tensor_selected(4);
M_rp = moment_tensor_selected(5);
M_tp = moment_tensor_selected(6);

% Seismic moment: M0 = sqrt(M_rr^2 + M_tt^2 + M_pp^2 + 2*(M_rt^2 + M_rp^2 + M_tp^2))
M0 = sqrt(M_rr^2 + M_tt^2 + M_pp^2 + 2*(M_rt^2 + M_rp^2 + M_tp^2));
% Convert from dyne·cm to N·m (conversion factor: 1e-7)
M0_Nm = M0 * 1e-7;
disp(['Seismic Moment (M0) = ' num2str(M0_Nm) ' N·m'])

%% Converting Epicentral Distance to km and Displaying Station Distances
distance_km = deg2km(distance);
for i = 1:length(station)
    fprintf('%s %.2f km\n', station{i}, distance_km(i));
end

%% Correct focal mechanism
% Assume moment_tensor_selected is a 6x1 vector:
% [mrr; mtt; mpp; mrt; mrp; mtp]
M_vec = moment_tensor_selected;  % e.g. 1e+25 * [7.3762; 9.7360; 4.9767; 7.5904; -5.7891; -5.7326]

% Extract components:
mrr = M_vec(1);
mtt = M_vec(2);
mpp = M_vec(3);
mrt = M_vec(4);
mrp = M_vec(5);
mtp = M_vec(6);

% Compute trace and isotropic component:
trace_val = mrr + mtt + mpp;
iso = trace_val / 3;

% Subtract isotropic part from the diagonal elements:
mrr_d = mrr - iso;
mtt_d = mtt - iso;
mpp_d = mpp - iso;

% Create the deviatoric moment tensor vector:
moment_tensor_deviatoric = [mrr_d, mtt_d, mpp_d, mrt, mrp, mtp];

% Plot the combined focal mechanism (beachball)
figure;
plotmt(1, 1, moment_tensor_deviatoric);
colormap([1 1 1; 1 0 0]);  % Typically red indicates compressional quadrants
title('The Focal Mechanism');





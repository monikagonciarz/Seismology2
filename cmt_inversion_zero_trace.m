% Matlab script with a least-squares algorithm to carry out a CMT source inversion.

% GEOL3048 Seismology II Practical week 9, Michael Frietsch and Ana Ferreira, March 2015
close all
clear all

% Load files with six sensitivity kernels, synthetic seismograms and real data in the data/ directory
folder_name='data';
cd(folder_name)
filenames={'mrr', 'mtt', 'mpp', 'mrt', 'mrp', 'mtp', 'real', 'synth'};
for n=1:length(filenames)
    load(filenames{n})
end
cd ..

% Inversion
% Decide input data (synthetics or real data) 
 inp_waves=real;
% inp_waves=synth;

% Build G matrix with the sensitivity kernels (sensitivity in u wrt each parameter), connect given data with framework for inversions
G = [mrr mtt mpp mrt mrp mtp];

% Impose the trace of the seismic moment tensor to be zero:
one = max(inp_waves)*1;
G = [G; [one one one 0 0 0]];
inp_waves = [inp_waves; 0];

% Carry out the source inversion using the input data
moment_tensor=inv(transpose(G)*G)*transpose(G)*inp_waves*10^26;

best_fit=[mrr mtt mpp mrt mrp mtp]*moment_tensor/10^26;

% Find the best double-couple fault geometry:
mainplane=mt2sdr(transpose(moment_tensor))
auxplane=auxplane(mainplane)





%% Plotting station distribution
% Load station file
fileID = fopen([folder_name '/stations']);
C = textscan(fileID, '%s %s %f %f');
station = C{1};
channel = C{2};
azimuth = C{3};
distance = C{4};

% Polar plot of the source and stations
figure(1)
h = polar(deg2rad(azimuth(1:1:end)),distance(1:1:end),'^');
hold on
set(h,'markersize',9,'MarkerFaceColor','b')
xlabel('Distance from Epicenter in Degree')
ylabel('Azimuth from North in Degree')
title('Station Distribution')

view(90,-90)
h = polar(0,0,'pr');
set(h,'markersize',15,'MarkerFaceColor','r')

offset = 3;
for i = 1:length(azimuth)
    angle_rad = deg2rad(azimuth(i));       % convert to radians
    [x, y] = pol2cart(angle_rad, distance(i) + offset); % convert polar to Cartesian
    % Create an outline by plotting the text several times in white, each slightly offset to simulate a border.
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





%% Plotting the comparison between input data and inversions results
load([folder_name '/cutpoints'])

figure(2)
[az, az_index]=sort(azimuth); % Plot with ascending azimuth
for m=1:length(station)
    k=az_index(m);% Plot with ascending azimuth
    subplot(7,3,m)
    window = cutpoints(k):cutpoints(k+1);
    time = 0:cutpoints(k+1)-cutpoints(k);

    % Calculate trace duration
    trace_duration = cutpoints(k+1) - cutpoints(k);
    fprintf('Trace duration for station %s: %.2f seconds\n', station{k}, trace_duration);

    plot(time, inp_waves(window),'k', time, best_fit(window),'r')
    xlim([0 800])
    ylim([-abs(max(inp_waves)) abs(max(inp_waves))]) % Same scale for all traces
    title([station{k} ' ' channel{k}(3) ' dist ' num2str(distance(k)) ' az ' num2str(azimuth(k))])
end

% Add a main title for the entire figure
sgtitle('Waveform Comparison for All Stations')
legend('Synthetic Data', 'Inversion')





%% Plot the focal mechanism:
figure(3)
plotmt(1,1,transpose(moment_tensor))
colormap([1 1 1; 1 0 0]) % Red for compressional, white for dilatational
title('The Focal Mechanism')





%% Plotting an example of moment tensor component kernel 
figure(4)
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




%% Plotting sensitivity with respect to component of moment tensor
figure(5)
[az, az_index]=sort(azimuth); % Plot with ascending azimuth
for m=1:length(station)
    k=az_index(m);% Plot with ascending azimuth
    subplot(7,3,m)
    window = cutpoints(k):cutpoints(k+1);
    time = 0:cutpoints(k+1)-cutpoints(k);
    plot(time, mtp(window),'b')
    xlim([0 800])
    ylim([-abs(max(mtp)) abs(max(mtp))]) % Same scale for all traces
    title([station{k} ' ' channel{k}(3) ' dist ' num2str(distance(k)) ' az ' num2str(azimuth(k))])
end
sgtitle('Sensitivity with respect to component of moment tensor')




%% Calculating the trace of the moment tensor
% trace = mrr + mtt + mpp
trace = moment_tensor(1,:) + moment_tensor(2,:) + moment_tensor(3,:)





%% An eigenvalue-eigenvector decomposition of the seismic moment tensor
M = [moment_tensor(1,:) moment_tensor(4,:) moment_tensor(5,:); moment_tensor(4,:) moment_tensor(2,:) moment_tensor(6,:); moment_tensor(5,:) moment_tensor(6,:) moment_tensor(3,:)];
eig(M)





%% Calculating the percentage deviation from the pure double couple component
[V,D]=eig(M);
epsilon = 100 * min(abs(diag(D))) / max(abs(diag(D)))





%% Calculating the seismic moment
M_rr = moment_tensor(1,:);
M_tt = moment_tensor(2,:);
M_pp = moment_tensor(3,:);
M_rt = moment_tensor(4,:);
M_rp = moment_tensor(5,:);
M_tp = moment_tensor(6,:);

% Calculate the seismic moment
M0 = sqrt(M_rr.^2 + M_tt.^2 + M_pp.^2 + 2*(M_rt.^2 + M_rp.^2 + M_tp.^2)); 

% Convert from dyne·cm to N·m (seismic moment in dyne·cm)
M0_Nm = M0 * 10^(-7);  % Conversion factor to N·m

% Output the seismic moment (in N·m)
disp(['Seismic Moment (M0) = ', num2str(M0_Nm), ' N·m'])





%% Converting epicentral distance to km
distance_km = deg2km(distance);

% Display the station names and their corresponding distances in kilometers
for i = 1:length(station)
    fprintf('%s %.2f km\n', station{i}, distance_km(i));
end


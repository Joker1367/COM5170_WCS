clc; clear; close all;
% Parameters
rng(1000);
N = 16;              % Antenna number
d = 0.5;             % Antenna spacing
L = 5;               % Singal number
angles = rand(1, L) * pi - pi/2;  % incidence angles [-π/2, π/2]

pic_index = 1;

%% Q1
% Steering Vectors
steering_vectors = zeros(N, L);
for l = 1 : L
    theta = angles(l);
    steering_vectors(:, l) = exp(-1j * 2 * pi * d * (0:N-1)' * cos(theta));
end

disp('Steering Vectors:');
for l = 1 : L
    fprintf('Signal %d: Angle = %.2f degree\n', l, 180 * angles(l) / pi);
    disp(steering_vectors(:, l).');
end

%% Q2
azimuth_range = linspace(-pi, pi, 10000);  % angle [-π, π]
radiation_patterns = zeros(L, length(azimuth_range));

% radiation pattern
for l = 1 : L
    w = steering_vectors(:, l);  % beamforming vector
    for t = 1:length(azimuth_range)
        phi = azimuth_range(t);  
        a_phi = exp(-1j * 2 * pi * d * (0:N-1)' * cos(phi));  
        radiation_patterns(l, t) = abs(w' * a_phi);      
    end
end

% plot radiation pattern
for l = 1 : L
    figure(pic_index);
    polarplot(azimuth_range, radiation_patterns(l, :));
    title(sprintf('Radiation Pattern for Signal %d (phi = %.2f degree)', l, 180 * angles(l) / pi));
    pic_index = pic_index + 1;
end

%% Q3
% Parameters
P = 1;  % Power of each signal, assumed to be 1
SIR_values = zeros(1, L);  % To store SIR for each signal

for target_signal = 1:L
    % Beamforming vector for the target signal
    w = steering_vectors(:, target_signal); 
    
    % Desired signal power
    desired_signal_power = P * abs(w' * steering_vectors(:, target_signal))^2;
    
    % Interference signal power
    interference_power = 0;
    for interference_signal = 1:L
        if interference_signal ~= target_signal
            interference_power = interference_power + ...
                P * abs(w' * steering_vectors(:, interference_signal))^2;
        end
    end
    
    % Calculate SIR
    SIR_values(target_signal) = desired_signal_power / interference_power;
end

% SIR in decibels
SIR_dB = 10 * log10(SIR_values);

% Display results
disp(sprintf('SIR (N = %d) values for each signal:', N));
disp(SIR_dB);

%% Q4
N = 64;
steering_vectors = zeros(N, L);
for l = 1 : L
    theta = angles(l);
    steering_vectors(:, l) = exp(-1j * 2 * pi * d * (0:N-1)' * cos(theta));
end

P = 1;  % Power of each signal, assumed to be 1
SIR_values = zeros(1, L);  % To store SIR for each signal

for target_signal = 1:L
    % Beamforming vector for the target signal
    w = steering_vectors(:, target_signal); 
    
    % Desired signal power
    desired_signal_power = P * abs(w' * steering_vectors(:, target_signal))^2;
    
    % Interference signal power
    interference_power = 0;
    for interference_signal = 1:L
        if interference_signal ~= target_signal
            interference_power = interference_power + ...
                P * abs(w' * steering_vectors(:, interference_signal))^2;
        end
    end
    
    % Calculate SIR
    SIR_values(target_signal) = desired_signal_power / interference_power;
end

% SIR in decibels
SIR_dB = 10 * log10(SIR_values);

% Display results
disp(sprintf('SIR (N = %d) values for each signal:', N));
disp(SIR_dB);

azimuth_range = linspace(-pi, pi, 10000);  % angle [-π, π]
radiation_patterns = zeros(L, length(azimuth_range));

% radiation pattern
for l = 1 : L
    w = steering_vectors(:, l);  % beamforming vector
    for t = 1:length(azimuth_range)
        phi = azimuth_range(t);  
        a_phi = exp(-1j * 2 * pi * d * (0:N-1)' * cos(phi));  
        radiation_patterns(l, t) = abs(w' * a_phi);      
    end
end

% plot radiation pattern
for l = 1 : L
    figure(pic_index);
    polarplot(azimuth_range, radiation_patterns(l, :));
    title(sprintf('Radiation Pattern for Signal %d (phi = %.2f degree and Nr = 64)', l, 180 * angles(l) / pi));
    pic_index = pic_index + 1;
end
clc;
clear;
close all;



% Basic parameters
K = 1;                          % Rician K-factor
s = 1;                          % direct path
L_vals = [1 2 3 4];             % Number of diversity branches
Es_N0_dB = [1, 3, 5, 7, 9];     % Symbol energy-to-noise ratio (in dB)
Es_N0 = 10.^(Es_N0_dB / 10);    % Linear scale
numBits = 1e6;                  % Total number of bits for simulation
Pb_target = 1e-4;               % Target BER

signal_real = 2 * randi([0, 1], 1, numBits) - 1;
signal_imag = 2 * randi([0, 1], 1, numBits) - 1;
signal = (1/sqrt(2)) * (signal_real + 1i * signal_imag);

% Initialize BER
ber_SC = zeros(length(L_vals), length(Es_N0_dB));
ber_MRC = zeros(length(L_vals), length(Es_N0_dB));
ber_EGC = zeros(length(L_vals), length(Es_N0_dB));
ber_DC = zeros(length(L_vals), length(Es_N0_dB));

% Main loop: Iterate through different numbers of diversity branches
for l_idx = 1:length(L_vals)
    L = L_vals(l_idx); % Number of diversity branches
    
    for snr_idx = 1:length(Es_N0_dB)
        snr = Es_N0(snr_idx); % Linear SNR

        % Multiple branch Rician fading
        h = sqrt(K / (K + 1)) * s + sqrt(1 / (K + 1)) * (1/sqrt(2)) * (randn(L, numBits) + 1i*randn(L, numBits)); % Rayleigh channel for each branch
        noise = sqrt(1/(2*snr)) * (1/sqrt(2)) * (randn(L, numBits) + 1i*randn(L, numBits)); % AWGN
        
        % Process different diversity combining techniques
        r = zeros(L, numBits);
        for i = 1:L
            r(i, :) = h(i, :) .* signal + noise(i, :);
        end
        
        % (a) Selective Combining (SC)
        [max_gain, idx] = max(abs(r), [], 1); % Select the branch with maximum gain
        demod_sc = r(idx + (0:size(r, 2)-1)*size(r, 1)) ./ h(idx + (0:size(r, 2)-1)*size(r, 1));
        bits_rx_sc_real = sign(real(demod_sc));
        bits_rx_sc_imag = sign(imag(demod_sc));
        ber_SC(l_idx, snr_idx) = (sum(bits_rx_sc_real~=signal_real) + sum(bits_rx_sc_imag~=signal_imag)) / (2 * numBits);
        
        % (b) Maximal Ratio Combining (MRC)
        r_mrc = sum(conj(h) .* r, 1); % MRC combining
        demod_mrc = r_mrc; % Demodulation
        bits_rx_mrc_real = sign(real(demod_mrc));
        bits_rx_mrc_imag = sign(imag(demod_mrc));
        ber_MRC(l_idx, snr_idx) = (sum(bits_rx_mrc_real~=signal_real) + sum(bits_rx_mrc_imag~=signal_imag)) / (2 * numBits);
        
        % (c) Equal Gain Combining (EGC)
        r_egc = sum(r .* exp(-1i * angle(h)), 1); % EGC combining
        demod_egc = r_egc; % Demodulation
        bits_rx_egc_real = sign(real(demod_egc));
        bits_rx_egc_imag = sign(imag(demod_egc));
        ber_EGC(l_idx, snr_idx) = (sum(bits_rx_egc_real~=signal_real) + sum(bits_rx_egc_imag~=signal_imag)) / (2 * numBits);
        
        % (d) Direct Combining (DC)
        r_dc = sum(r, 1); % Direct combining
        r_dc = r_dc .* exp(-1i * angle(sum(h, 1))); % Phase correction
        demod_dc = r_dc; % Demodulation
        bits_rx_dc_real = sign(real(demod_dc));
        bits_rx_dc_imag = sign(imag(demod_dc));
        ber_DC(l_idx, snr_idx) = (sum(bits_rx_dc_real~=signal_real) + sum(bits_rx_dc_imag~=signal_imag)) / (2 * numBits);
    end
end

% Plotting: Four figures for the four diversity combining techniques
markers = {'o-', 's-', 'd-', '^-'}; % Different markers

% (a) Selective Combining
figure;
for l_idx = 1:length(L_vals)
    semilogy(Es_N0_dB, ber_SC(l_idx, :), markers{l_idx}, 'DisplayName', sprintf('L=%d', L_vals(l_idx)));
    hold on;
end
grid on;
xlabel('E_s/N_0 (dB)');
ylabel('Bit Error Rate (BER)');
ylim([1e-6, 1]); % Set y-axis limits
legend show;
title('Selective Combining (SC, Rician)');

% (b) Maximal Ratio Combining
figure;
for l_idx = 1:length(L_vals)
    semilogy(Es_N0_dB, ber_MRC(l_idx, :), markers{l_idx}, 'DisplayName', sprintf('L=%d', L_vals(l_idx)));
    hold on;
end
grid on;
xlabel('E_s/N_0 (dB)');
ylabel('Bit Error Rate (BER)');
ylim([1e-6, 1]); % Set y-axis limits
legend show;
title('Maximal Ratio Combining (MRC, Rician)');

% (c) Equal Gain Combining
figure;
for l_idx = 1:length(L_vals)
    semilogy(Es_N0_dB, ber_EGC(l_idx, :), markers{l_idx}, 'DisplayName', sprintf('L=%d', L_vals(l_idx)));
    hold on;
end
grid on;
xlabel('E_s/N_0 (dB)');
ylabel('Bit Error Rate (BER)');
ylim([1e-6, 1]); % Set y-axis limits
legend show;
title('Equal Gain Combining (EGC, Rician)');

% (d) Direct Combining
figure;
for l_idx = 1:length(L_vals)
    semilogy(Es_N0_dB, ber_DC(l_idx, :), markers{l_idx}, 'DisplayName', sprintf('L=%d', L_vals(l_idx)));
    hold on;
end
grid on;
xlabel('E_s/N_0 (dB)');
ylabel('Bit Error Rate (BER)');
ylim([1e-6, 1]); % Set y-axis limits
legend show;
title('Direct Combining (DC, Rician)');


clc;
clear all;
close all;

%% Parameters
num_samples = 1000;
sampling_rate = 1000; % in samples per second
symbols_per_frame = 10; % Number of symbols per frame
num_subcarriers = symbols_per_frame; % Number of subcarriers
subcarrier_spacing = sampling_rate / symbols_per_frame; % Spacing between subcarriers
%carrier_freq = 10 * subcarrier_spacing; % Carrier frequency (in Hz)
M=16;

%% QAM constellation size

random_symbols = randi([0 M-1], 1, num_samples);
qam_symbols = qammod(random_symbols, M);

%% Serial to parallel mapping
num_frames = ceil(num_samples / symbols_per_frame);
padded_data = [qam_symbols, zeros(1, num_frames * symbols_per_frame - num_samples)]; % Padding if needed
parallel_data = reshape(padded_data, symbols_per_frame, num_frames).'; % Reshape to matrix

% Plot the parallel data
figure;
plot(real(parallel_data(:)), imag(parallel_data(:)), 'x');
xlabel('I');
ylabel('Q');
title('QAM Symbols');
axis square;

%% Generate time vector
t_frame = linspace(0, 1/sampling_rate, symbols_per_frame); % Time for one frame
t = repelem(t_frame, 1, num_frames); % Replicate for all frames

%% Perform IFFT for each frame
ofdm_signals = zeros(num_frames, symbols_per_frame+4);
for i = 1:num_frames
    % Apply IFFT
    ofdm_signal_time_domain = ifft(parallel_data(i, :), num_subcarriers);
    % Add cyclic prefix (optional but common in OFDM)
    cp_length = 4; % Length of cyclic prefix
    ofdm_signal_cp = [ofdm_signal_time_domain(end-cp_length+1:end), ofdm_signal_time_domain(1:end)]; % Append cyclic prefix
    % Store the OFDM signal for this frame
    ofdm_signals(i, :) = ofdm_signal_cp;
end

% Plot the OFDM signals
figure;
imagesc(real(ofdm_signals));
xlabel('Sample Index');
ylabel('Frame Index');
title('OFDM Signals');
colorbar;

%% Parallel to Serial
% Convert parallel data to serial data
serial_data = reshape(ofdm_signals.', 1, []);

% Adjust the time vector to match the length of serial data
t_serial = linspace(0, (num_samples / sampling_rate), length(serial_data));

% Plot the serial data (only the real part)
figure;
plot(t_serial, real(serial_data));
xlabel('Time (s)');
ylabel('Amplitude');
title('Serial Data (Real Part)');
grid on;

%% Convert digital to analog (DA) - Interpolation
da_factor = 10; % Upsampling factor
t_da = linspace(0, max(t_serial), length(serial_data) * da_factor); % Time vector for DA
serial_data_da = interp1(t_serial, serial_data, t_da, 'linear'); % Upsample using linear interpolation

% Plot the digital signal before and after conversion to analog
figure;
subplot(2, 1, 1);
plot(t_serial, real(serial_data));
xlabel('Time (s)');
ylabel('Amplitude');
title('Digital Signal (Real Part)');
grid on;

subplot(2, 1, 2);
plot(t_da, real(serial_data_da));
xlabel('Time (s)');
ylabel('Amplitude');
title('Analog Signal (Real Part)');
grid on;

%% Modulate the analog signal (BPSK modulation)
carrier_freq = 1000; % Carrier frequency for modulation
analog_signal = real(serial_data_da) .* cos(2 * pi * carrier_freq * t_da); % BPSK modulation

% Plot the modulated signal
figure;
plot(t_da, analog_signal);
xlabel('Time (s)');
ylabel('Amplitude');
title('BPSK Modulated Signal');
grid on;

%% Add noise to the modulated signal
SNR_dB = 1000; % Signal-to-Noise Ratio (in dB)
analog_signal_noisy = awgn(analog_signal, SNR_dB, 'measured');

% Plot the noisy signal
figure;
plot(t_da, analog_signal_noisy);
xlabel('Time (s)');
ylabel('Amplitude');
title('Noisy Analog Signal (BPSK Modulated)');
grid on;

%% BPSK Demodulation
demodulated_signal = real(analog_signal_noisy) .* cos(2 * pi * carrier_freq * t_da);

% Plot the demodulated signal
figure;
plot(t_da, demodulated_signal);
xlabel('Time (s)');
ylabel('Amplitude');
title('Demodulated Signal (Before Filtering)');
grid on;

%% Low-pass filter to remove high-frequency noise
% Correct the cutoff frequency calculation
cutoff_freq = 0.2; % Cutoff frequency as a fraction of the Nyquist frequency

% Design the Butterworth filter
[b, a] = butter(6, cutoff_freq, 'low');
 % Design the Butterworth filter
demodulated_signal_filtered = filter(b, a, demodulated_signal); % Apply the filter

figure;
plot(t_da, demodulated_signal_filtered);
xlabel('Time (s)');
ylabel('Amplitude');
title('Demodulated Signal (After Filtering)');
grid on;

%% Analog to digital conversion
num_received_samples = length(demodulated_signal_filtered); % Adjusted number of samples after filtering
t_received = linspace(0, (num_received_samples / sampling_rate), num_received_samples); % Time vector for received signal

% Downsample to match the original sampling rate
downsample_factor = length(serial_data) / num_received_samples;
serial_data_digital_received = interp1(t_received, demodulated_signal_filtered, t_serial, 'linear');

% Plot the digital signal before and after conversion from analog
figure;
plot(t_serial, real(serial_data));
hold on;
plot(t_serial, real(serial_data_digital_received), 'r--');
hold off;
xlabel('Time (s)');
ylabel('Amplitude');
title('Digital Signal (Before and After Conversion)');
legend('Original Digital Signal', 'Recovered Digital Signal');
grid on;

%% Convert serial data to parallel data
received_ofdm_signals = reshape(serial_data_digital_received, symbols_per_frame + cp_length, []);

%% Remove cyclic prefix and perform FFT for each frame
received_data_freq_domain = zeros(num_frames, symbols_per_frame);
for i = 1:num_frames
    % Remove cyclic prefix
    received_ofdm_signal_no_prefix = received_ofdm_signals(cp_length+1:end, i);
    % Perform FFT
    received_data_freq_domain(i, :) = fft(received_ofdm_signal_no_prefix, num_subcarriers);
end

% Plot the received data in frequency domain
figure;
subplot(3, 1, 1);
imagesc(abs(received_data_freq_domain));
xlabel('Subcarrier Index');
ylabel('Frame Index');
title('Received Data (Frequency Domain)');
colorbar;

% Plot the phase of the received data in frequency domain
subplot(3, 1, 2);
imagesc(angle(received_data_freq_domain));
xlabel('Subcarrier Index');
ylabel('Frame Index');
title('Phase of Received Data (Frequency Domain)');
colorbar;

% Plot the magnitude of the received data in frequency domain
subplot(3, 1, 3);
imagesc(log(abs(received_data_freq_domain) + 1)); % Adding 1 to avoid log(0)
xlabel('Subcarrier Index');
ylabel('Frame Index');
title('Log Magnitude of Received Data (Frequency Domain)');
colorbar;


%% Convert parallel data from FFT output to serial data
received_serial_data = reshape(received_data_freq_domain.', 1, []);

% Plot the received serial data
figure;
plot(real(received_serial_data));
xlabel('Sample Index');
ylabel('Amplitude');
title('Received Serial Data (Real Part)');
grid on;


%% QAM Demodulation
received_qam_symbols = qamdemod(received_serial_data, M);

%%Plot the received QAM symbols
figure;
plot(real(received_qam_symbols), imag(received_qam_symbols), 'x');
xlabel('I');
ylabel('Q');
title('Received QAM Symbols');
axis square;





%% BER vs SNR Curve

% Initialize arrays to store BER values
SNR_range = -10:2:20; % Range of SNR values to simulate
BER = zeros(size(SNR_range));

for i = 1:length(SNR_range)
    % Add noise to the modulated signal
    SNR_dB = SNR_range(i); % Current SNR value
    analog_signal_noisy = awgn(analog_signal, SNR_dB, 'measured');
    
    % BPSK Demodulation
    demodulated_signal = real(analog_signal_noisy) .* cos(2 * pi * carrier_freq * t_da);
    
    % Low-pass filter to remove high-frequency noise
    demodulated_signal_filtered = filter(b, a, demodulated_signal); % Apply the filter
    
    % Analog to digital conversion
    num_received_samples = length(demodulated_signal_filtered); % Adjusted number of samples after filtering
    t_received = linspace(0, (num_received_samples / sampling_rate), num_received_samples); % Time vector for received signal
    
    % Downsample to match the original sampling rate
    serial_data_digital_received = interp1(t_received, demodulated_signal_filtered, t_serial, 'linear');
    
    % Convert serial data to parallel data
    received_ofdm_signals = reshape(serial_data_digital_received, symbols_per_frame + cp_length, []);
    
    % Remove cyclic prefix and perform FFT for each frame
    received_data_freq_domain = zeros(num_frames, symbols_per_frame);
    for j = 1:num_frames
        % Remove cyclic prefix
        received_ofdm_signal_no_prefix = received_ofdm_signals(cp_length+1:end, j);
        % Perform FFT
        received_data_freq_domain(j, :) = fft(received_ofdm_signal_no_prefix, num_subcarriers);
    end
    
    % Convert parallel data from FFT output to serial data
    received_serial_data = reshape(received_data_freq_domain.', 1, []);
    
    % QAM Demodulation
    received_qam_symbols = qamdemod(received_serial_data, M);
    
    % Calculate BER
    BER(i) = sum(received_qam_symbols ~= random_symbols) / length(random_symbols);
end

% Plot BER vs SNR curve
figure;
semilogy(SNR_range, BER, 'bo-');
xlabel('SNR (dB)');
ylabel('Bit Error Rate (BER)');
title('BER vs SNR');
grid on;

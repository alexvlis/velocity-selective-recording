function out_sig = electroshift(signal, shifts, velocity)
% Applies pulse shift in time due to electrode spacing errors.
%   
% A Vlissidis

fs = 5e5; % Sampling Frequency
electro_length = 1e-3; % Electrode spacing

[len, numchannels] = size(signal);

% Copy the time series and first channel data verbatim
out_sig = zeros(len, numchannels-1);
out_sig(:, 1:2) = signal(:, 1:2);

% Shift each channel
for i = 3:numchannels
    % Calculate the time shift
    t = (shifts(i-2) - electro_length)/velocity;
    % Translate time to number of samples
    N = fs * t;
    out_sig(:, i) = circshift(signal(:, i), round(N));
end

% Plot the shifted signal
plot_mec_signal(out_sig, 1, numel(signal(:, 1)));

end
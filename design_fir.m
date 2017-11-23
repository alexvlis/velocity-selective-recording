function design_fir(Parameters, velocity, v2)
% Design FIR filters to perform VSR
%   
% A. Vlissidis

% Generate the input signal x
velocities = 10:1:120;
InitTime = 0.001;
EndTime = 0.004;
time = -InitTime:1/Parameters.SamplingFrequency:EndTime;
timelen = numel(time);
global_time = (0:(timelen * numel(velocities)) - 1) ...
                   / Parameters.SamplingFrequency;

data = zeros(numel(velocities), 3, timelen);
for i = 1:numel(velocities)
    % Get the velocity
    v = velocities(i);
    
    % Calculate the time range
    start = 1 + (i - 1) * timelen;
    tend = start + timelen - 1;
    
    % Insert the AP into the input sequence
    x(:, start:tend) = GetBiPolar(Parameters, v, time);
    data(i, :, :) = GetBiPolar(Parameters, v, time);
end
data = AgcSim(data);
y = zeros(size(x)); % Allocate the space for the output
    
% Initialize the coefficients for the filters
a = 1; % FIR
b = zeros(Parameters.channels, Parameters.AnnLength(1)); % Preallocate for speed
c = zeros(Parameters.channels, Parameters.AnnLength(1));

Plot the input
figure;
for channel = 1:Parameters.channels
    subplot(Parameters.channels, 1, channel);
    plot(global_time, x(channel, :))
    if channel == 1
       title('Input Action Potential') 
    end
end
xlabel('time (s)')

% Calculate the interchannel delay in number of samples
n = Parameters.ElectrodeSpacing * 2 * Parameters.SamplingFrequency/velocity;
n2 = Parameters.ElectrodeSpacing * 2 * Parameters.SamplingFrequency/v2;

figure;
hold on
for channel = 1:Parameters.channels
    % Apply centroid
    b(channel, :) = -2 * (1:Parameters.AnnLength(1))/Parameters.AnnLength(1) + 1;
    x(channel, :) = filter(b(channel, :), a, x(channel, :));
    
    if channel < Parameters.channels
        % Apply delay
        N = (Parameters.channels - channel) * round(n);
        N2 = (Parameters.channels - channel) * round(n2);
        c(channel, :) = 0;
        if mod(channel, 2) == 0
            c(channel, N2) = -1.803;
        else 
            c(channel, N) = 0.3;
        end
        y(channel, :) = filter(c(channel, :), a, x(channel, :));
    else
        y(channel, :) = 1.813 * x(channel, :); % scale last channel
    end
    % Plot the result
    subplot(Parameters.channels, 1, channel);
    plot(global_time, y(channel, :))
    if channel == 1
      title('FIR Outputs') 
    end
end
xlabel('time (s)')
legend('1', '2', '3')
hold off

% Sum the output of each channel
h = sum(y);
figure;
plot(global_time, h);
title('Summed Output');
xlabel('Time (s)');
ylabel('Amplitude');

peaks = zeros(1, numel(velocities));
for i = 1:numel(velocities)
    % Get the velocity
    h = sim_fir(b, c, Parameters, reshape(data(i, :, :), 3, 501));
    peaks(i) = max(abs(h));
end
% Normalise
peaks = peaks./max(peaks);

figure;
plot(velocities, peaks)
title('Velocity Spectrum')
xlabel('Velocity m/s')
ylabel('Magnitude of Response')
end
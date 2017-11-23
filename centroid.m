% Setup the FIR filter parameters
N = 50;
a = 1;
b = -2 * (1:N)./N + 1;

% Generate the input
InitTime = 0.001;
EndTime = 0.004;
time = -InitTime:1/Parameters.SamplingFrequency:EndTime;
x = GetBiPolar(Parameters, 50, time);

y = filter(b, a, x(1, :));

figure;
hold on
plot(time*1000, 8*x(1, :), '--');
plot(time*1000, y)
title('Centroid Filter Response')
xlabel('Time (ms)')
ylabel('Response')
legend('Input AP', 'Filter Response')
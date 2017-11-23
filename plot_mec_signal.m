function plot_mec_signal(data, start, stop)
%PLOT_RAW_DATA Summary of this function goes here
%   Detailed explanation goes here

[~, channels] = size(data);

figure;
for i = 2:channels
    subplot(channels-1, 1, i-1);
    plot(data(start:stop, 1), data(start:stop, i));
    if i == 2
       title('L5 Dermatome Stimulation'); 
    end
end
xlabel('Time (s)');

end


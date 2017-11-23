function sim_real_data(ANN, rawdata, start, stop)
% Test the model with real data
%
% A Vlissidis

matched = [10, 20, 30];
%matched = 5:5:30;

fs = 1e5;
fc = [200 5000];

figure;
for i = 2:6
    subplot(5, 1, i-1);
    plot(rawdata(start:5:stop, 1),rawdata(start:5:stop, i));
end
xlabel('Time (s)');

% Apply electro shifts
data = electroshift(rawdata(start:stop, :), [0.7208 1.2092 0.8901 1.1302] * 1e-3, 8);
data = data(1:5:end, :);
InputSequence = data(:, 2:end)';

% Filter the data to remove Red Noise
%medFilt = dsp.MedianFilter(2);
%InputSequence = medFilt(InputSequence);
% [b, a] = butter(3, fc/(fs/2), 'bandpass');
% [channels, ~] = size(InputSequence);
% for i = 1:channels
%    InputSequence(i, :) = filter(b, a, InputSequence(i, :)); 
% end

% Plot the filtered input channels
% figure;
% for i = 2:6
%     subplot(5, 1, i-1);
%     plot(rawdata(start:5:stop, 1),InputSequence(i-1, :));
% end
% xlabel('Time (s)');

NormInput = AgcSim(InputSequence);
InputCellArray = con2seq(NormInput);

figure;
NumAnns = numel(ANN);
buf = zeros(1, NumAnns);
for i = 1:NumAnns
   net = ANN(i).net;
   OutputCellArray = sim(net, InputCellArray);
   AnnOutput = cell2mat(OutputCellArray);
   buf(i) = max(AnnOutput);
   subplot(numel(ANN), 1, i);
   plot(rawdata(start:5:stop, 1), AnnOutput);
end
xlabel('Time (s)')

[~, idx] = max(buf);
matched(idx)

end


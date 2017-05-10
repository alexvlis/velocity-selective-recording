%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Demonstrate velocity resolution of delay and add method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plots corresponding velocity against the number of samples delayed

Fs = 100000;    % Sampling frequency
d = 0.003;      % Electrode spacing
maxsamples = 15;
v = zeros(maxsamples);

for samples = 1:maxsamples
    v(samples) = d/(samples/Fs);
end

figure;
stem(v, 'bo');
ax = gca;
ax.XTick = 0:maxsamples;
grid on;
xlim([0 maxsamples]);
xlabel(['Number of samples delayed (F_s = ' num2str(Fs/1000) 'kHz)']);
ylabel('Velocity (m/s)');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Display unipolar TMAPs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Parameters = struct (                 ...
   'Electrodes',              11    , ...
   'ElectrodeSpacing',      3e-3    , ...
   'SamplingFrequency',      1e5    , ...
   'ActionPotentialType',    'tmap2', ...
   'StartTestVelocity',       10    , ...
   'StepTestVelocity',         1    , ...
   'EndTestVelocity',        120    , ...
   'MatchRepeats',            10    , ...
   'NoiseLevel',               0    , ...
   'APType',                'UniPolar');

% Time sequence length for each AP
InitTime = 0.001;
EndTime = 0.004;
Time = -InitTime:1/Parameters.SamplingFrequency:EndTime;

% Sensor type
if Parameters.APType == 'UniPolar'
    GetData = @GetUniPolar;
    DataLines = Parameters.Electrodes;
elseif Parameters.APType == 'TriPolar'
    GetData = @GetTriPolar;
    DataLines = Parameters.Electrodes - 2;
end

% AP velocity
ActionPotentialVelocity = 20;
tmap_types = {'tmap1', 'tmap2'};

figure;
for i = 1:numel(tmap_types)
    Parameters.ActionPotentialType = tmap_types{i};
    UniPolarSignal = GetUniPolar(Parameters, ActionPotentialVelocity, Time);
    subplot(numel(tmap_types),1,i);
    plot(Time,UniPolarSignal(1,:));
    title(tmap_types(i));
end
xlabel('Time (s)');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Display tripolar APs and target pulses of different velocities
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Parameters = struct (                 ...
   'Electrodes',              11    , ...
   'ElectrodeSpacing',      3e-3    , ...
   'SamplingFrequency',      1e6    , ...
   'ActionPotentialType',    'tmap2', ...
   'StartTestVelocity',       20    , ...
   'StepTestVelocity',        90    , ...
   'EndTestVelocity',        100    , ...
   'MatchRepeats',            10    , ...
   'NoiseLevel',               0    );

InitTime = 0.001;
EndTime = 0.004;
Time = -InitTime:1/Parameters.SamplingFrequency:EndTime;

Velocities = Parameters.StartTestVelocity:...
             Parameters.StepTestVelocity:...
             Parameters.EndTestVelocity;
NumVelocities = numel(Velocities);
Data = zeros(NumVelocities , DataLines , numel(Time));

for VelocityIndex = 1:NumVelocities
    Velocity = Velocities(VelocityIndex);
    Data(VelocityIndex,:,:) = GetData(Parameters, Velocity, Time);
end

TripolarNorm = AgcSim(Data);
FinalChannel = TripolarNorm(:, DataLines,:);
TargetPulse = zeros(size(FinalChannel));
for i = 1:NumVelocities
    TargetPulse(i,1,(FinalChannel(i,1,:) >= (max(FinalChannel(i,1,:)/sqrt(2))))) = 1;
end

figure;
for n = 1:NumVelocities
    Velocity = Velocities(n);
    subplot(NumVelocities,1,n);
    hold on;
    plot(Time, reshape(TripolarNorm(n,9,:), [1 numel(Time)]));
    plot(Time, reshape(TargetPulse(n,1,:), [1 numel(Time)]));
    title([num2str(Velocity) ' m/s']);
    grid on;
    ylim([-1 1]);
end
xlabel('Time (s)');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Delay and sum method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Parameters = struct (                 ...
   'Electrodes',              11    , ...
   'ElectrodeSpacing',      3e-3    , ...
   'SamplingFrequency',      1e5    , ...
   'ActionPotentialType',    'tmap2', ...
   'StartTestVelocity',       10    , ...
   'StepTestVelocity',       0.1    , ...
   'EndTestVelocity',        160    , ...
   'NoiseLevel',               0    );

MatchedVelocity = [20, 50, 90];

% Set the initial velocity set including the repeated matched set
Velocities = Parameters.StartTestVelocity:...
             Parameters.StepTestVelocity:...
             Parameters.EndTestVelocity;
NumVelocities = numel(Velocities);

% Set up an output array
NumMatchedVelocities = numel(MatchedVelocity);
PeakResponse = zeros(NumVelocities,NumMatchedVelocities);

% Set up the time sequence for the signals at each velocity
InitTime = 0.001;
EndTime = 0.004;
Time = -InitTime:1/Parameters.SamplingFrequency:EndTime;
TimeLength = numel(Time);

% Set up an array for the input tripoles so that we can do an overall AGC
Data = zeros(NumVelocities, DataLines, TimeLength);
           
% Run a test at each of the selected velocities
for VelocityIndex = 1:NumVelocities     
    
    % Retrieve the selected velocity
    Velocity = Velocities(VelocityIndex);
    
    % Get a data set of the approriate velocity
    Data(VelocityIndex,:,:) = GetData(Parameters, Velocity, Time);
end

% Apply the AGC
TripolarNorm = AgcSim(Data);

% Apply each delay-sum system to an AP at each test velocity
% and find the peak response.

DelaySumOutput = zeros(NumVelocities, 1, TimeLength);

% Apply data of the selected velocity
for MatchedVelInd=1:NumMatchedVelocities
    
    % Calculate number of samples for delay compensation
    InterElectrodeDelay = Parameters.ElectrodeSpacing/...
                          MatchedVelocity(MatchedVelInd);
    for VelocityIndex = 1:NumVelocities
        
        % Delay and sum
        DelaySumOutput(VelocityIndex,1,:) = TripolarNorm(VelocityIndex,DataLines,:);
        for n = 1:DataLines-1
            Delay = round((DataLines-n)*...
                          InterElectrodeDelay*...
                          Parameters.SamplingFrequency);
            DelayedSignal = circshift(TripolarNorm(VelocityIndex,n,:), Delay, 3);
            DelaySumOutput(VelocityIndex,:,:) = DelaySumOutput(VelocityIndex,:,:) + DelayedSignal;
        end
        
        % Get the peak reponse
        PeakResponse(VelocityIndex,MatchedVelInd) = max(max(DelaySumOutput(VelocityIndex,:,:)));
    end
end

% Use the peak responses to plot an IVS

% Plot the peaks
figure;
title('Normalised peak responses of delay-and-sum');
xlabel('Velocity (m/s)');
ylabel('Magnitude of response');
ylim([0 1.2]);
hold all;

% plot each of the selected matched velocities
for MatchedVelInd=1:NumMatchedVelocities
    
    % Normalise IVS
    PeakResponse(:,MatchedVelInd) = PeakResponse(:,MatchedVelInd)/...
                                           max(PeakResponse(:,MatchedVelInd));
    
    % Calculate velocity selectivity
    v3 = Velocities(PeakResponse(:,MatchedVelInd) >= (...
              PeakResponse(Velocities == MatchedVelocity(MatchedVelInd),MatchedVelInd)...
              /sqrt(2)));
    BW = max(v3) - min(v3);
    Qv = MatchedVelocity(MatchedVelInd) / BW;
    
    % Plot the peak reponse
    plot (Velocities,PeakResponse(:,MatchedVelInd),...
          'DisplayName',[num2str(MatchedVelocity(MatchedVelInd)) ' m/s']);
    
    text(MatchedVelocity(MatchedVelInd)-5,...
         max(PeakResponse(:,MatchedVelInd))+0.05, ['Q_v = ' num2str(Qv)]);    
end

% Show the legend
legend('show');

% % Delay and sum with BPF inputs
% Bpf = fir2(120,[0 0.5 0.5 0.55 0.55 1],[0 0 1 1 0 0]);
% for n = 1:DataLines
%     Outbpf(n,:) = filter(Bpf,1,TripolarNorm(n,:));
% end
% DelaySumOutput = Outbpf(9,:);
% InterElectrodeDelay = Parameters.ElectrodeSpacing/MatchedVelocity;
% for n = 1:Parameters.Electrodes-3
%     Delay = round((DataLines-n)*...
%                   InterElectrodeDelay*...
%                   Parameters.SamplingFrequency);
%     DelayedSignal = circshift(Outbpf(n,1:end), Delay, 2);
%     DelaySumOutput = DelaySumOutput + DelayedSignal;
% end

% figure;
% subplot(3,1,1);
% plot(SequenceTime, TripolarNorm(1,:));
% subplot(3,1,2);
% plot(SequenceTime, TripolarNorm(end,:));
% % plot(SequenceTime, Outbpf(i,:));
% subplot(3,1,3);
% plot(SequenceTime, DelaySumOutput);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compare tripolar input to pulse target output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Parameters = struct (                 ...
   'Electrodes',              11    , ...
   'ElectrodeSpacing',      3e-3    , ...
   'SamplingFrequency',      1e5    , ...
   'ActionPotentialType',    'tmap2', ...
   'StartTestVelocity',       20    , ...
   'StepTestVelocity',        20    , ...
   'EndTestVelocity',        100    , ...
   'MatchRepeats',             1    , ...
   'NoiseLevel',           2e-21    , ...
   'APType',              'UniPolar');

MatchedVelocity = 20;

% Set the initial velocity set including the repeated matched set
Velocities = cat(2,Parameters.StartTestVelocity: ...
    Parameters.StepTestVelocity:Parameters.EndTestVelocity, ...
    MatchedVelocity*ones(1,Parameters.MatchRepeats));
NumVelocities = numel(Velocities);

% Set up the time sequence for the signals at each velocity
InitTime = 0.001;
EndTime = 0.004;
Time = -InitTime:1/Parameters.SamplingFrequency:EndTime;
TimeLength = numel(Time);
SequenceTime = (0:(TimeLength * NumVelocities) - 1) ...
               / Parameters.SamplingFrequency;

% Set up the arrays to hold the input and output data
TripolarInputSequence = zeros(DataLines, ...
                              NumVelocities * TimeLength);
TargetSequence = zeros(1,NumVelocities*TimeLength);

% Create input signals
% Create an input sequence for the ANN
for VelocityIndex = 1:NumVelocities
    
    % Retrieve the selected velocity
    Velocity = Velocities(VelocityIndex);
    
    % Get a data set of the approriate velocity
    Data = GetData(Parameters, Velocity, Time);
    
    % Insert the data into the input sequence
    InsertStart = 1 + (VelocityIndex - 1) * TimeLength;
    InsertEnd = InsertStart + TimeLength - 1;
    TripolarInputSequence (:,InsertStart:InsertEnd) = Data;
end

% Normalise the data into the range 1 to -1 for each channel separately
TripolarNorm = AgcSim(TripolarInputSequence);

% Create an output target signal
% Add a pulse for each velocity that is the same as the matched
% velocity
VelocityIndices = 1:NumVelocities;
for VelocityIndex = VelocityIndices(Velocities == MatchedVelocity)
    
    % Get the insertion point for the data in the input sequence
    InsertStart = 1 + (VelocityIndex - 1) * TimeLength;
    InsertEnd = InsertStart + TimeLength - 1;
    
    % Get a data set of the approriate velocity with no noise
    NoiseLessParameters = Parameters;
    NoiseLessParameters.NoiseLevel = 0;
    Data = GetData(NoiseLessParameters, MatchedVelocity, Time);
    
    % Get the final channel
    FinalChannel = Data(DataLines,:);
    
    % Normalise it
    FinalChannel = AgcSim(FinalChannel);
    
    % Create a pulse at the 3dB points
    TargetPulse = zeros(size(FinalChannel));
    TargetPulse(FinalChannel >= (max(FinalChannel)/sqrt(2))) = 1;
    
    % Apply it to the target sequence
    TargetSequence(InsertStart:InsertEnd) = TargetPulse;
end

figure;
for i = 1:9
    subplot(10,1,i);
    plot(SequenceTime, TripolarNorm(i,:));
    xlim([0 0.03]);
end
subplot(10,1,10);
plot(SequenceTime, TargetSequence);
ylim([-1 1]);
xlim([0 0.03]);
xlabel('Time (s)');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ANN time response
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Parameters.StartTestVelocity = 10;
Parameters.StepTestVelocity = 10;
Parameters.EndTestVelocity = 120;

figure;
for ANNindex = 1:3
    ANN(ANNindex).net.outputConnect = [1 0];
    ANN(ANNindex).net.layers{1}.transferFcn = 'purelin';
    % Set the initial velocity set including the repeated matched set
    Velocities = Parameters.StartTestVelocity:...
                 Parameters.StepTestVelocity:...
                 Parameters.EndTestVelocity;
    NumVelocities = numel(Velocities);
    
    % Set up the time sequence for the signals at each velocity
    InitTime = 0.001;
    EndTime = 0.004;
    Time = -InitTime:1/Parameters.SamplingFrequency:EndTime;
    TimeLength = numel(Time);
    SequenceTime = (0:(TimeLength * NumVelocities) - 1) ...
        / Parameters.SamplingFrequency;
    
    % Run a test at each of the selected velocities
    for VelocityIndex = 1:NumVelocities
        
        % Retrieve the selected velocity
        Velocity = Velocities(VelocityIndex);
        
        % Get a data set of the approriate velocity
        Data = GetData(Parameters, Velocity, Time);
        
        % Insert the data into the input sequence
        InsertStart = 1 + (VelocityIndex - 1) * TimeLength;
        InsertEnd = InsertStart + TimeLength - 1;
        TripolarInputSequence (:,InsertStart:InsertEnd) = Data;
    end
    
    % Apply the AGC
    TripolarNorm = AgcSim(TripolarInputSequence);
    
    % Convert data format to be input to the ANN from a R-by-TS matrix
    % to a 1-by-TS cell array of R-by-1 vectors
    TripolarCellArray = con2seq(TripolarNorm);
    
    % Run the ANN
    AnnOutputCellArray = sim(ANN(ANNindex).net,TripolarCellArray);
    
    % Convert ANN output back to a conventional array
    AnnOutput = cell2mat(AnnOutputCellArray);
    AnnOutput = AgcSim(AnnOutput);
    
    subplot(4,1,1);
    plot(SequenceTime, TripolarNorm(Parameters.Electrodes-2,:));
    xlim([0 NumVelocities*0.005]);
    ylabel('Input');
    subplot(4,1,ANNindex+1);
    plot(SequenceTime, AnnOutput);
    ylabel([num2str(AnnVelocities(ANNindex)) ' m/s']);
    ylim([-1 1]);
    xlim([0 NumVelocities*0.005]);
end
xlabel('Time (s)');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute frequency spectrum of tripolar AP signals
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
InitTime = 0.001;
EndTime = 0.004;
Time = -InitTime:1/Parameters.SamplingFrequency:EndTime;

Velocities = [20, 50, 90];
NumVelocities = numel(Velocities);
Data = zeros(NumVelocities , DataLines , numel(Time));

for VelocityIndex = 1:NumVelocities
    Velocity = Velocities(VelocityIndex);
    Data(VelocityIndex,:,:) = GetData(Parameters, Velocity, Time);
end

TripolarNorm = AgcSim(Data);

figure;
for n = 1:NumVelocities
    Velocity = Velocities(n);
    Signal = reshape(TripolarNorm(n,end,:), [1 numel(Time)]);
    L = 2*floor(numel(Signal)/2);
    Y = fft(Signal);
    P2 = abs(Y/L);
    P1 = P2(1:L/2+1);
    P1(1:end) = 2*P1(1:end);
    Frequency = Parameters.SamplingFrequency*(0:(L/2))/L;
    %figure;
    subplot(1,2,1);
    plot(Time, Signal);
    title(strcat(Parameters.APType, ' AP time signal'));
    xlabel('Time (s)');
    ylabel('Amplitude');
    hold on;
    subplot(1,2,2);
    plot(Frequency,P1);
    title('Single-Sided Amplitude Spectrum');
    xlabel('f (Hz)');
    ylabel('|P1(f)|');
    hold on;
    labels{n} = [num2str(Velocity) ' m/s'];
end
legend(labels);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot weights and compute frequency responses
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ANNindex = 3;

% Get weights
InputWeights = zeros(DataLines,120);
for Channel = 1:(DataLines)
    i = 0;
    for n = Channel:...
            (DataLines):...
            (numel(ANN(ANNindex).net.IW{1}) - (DataLines) + Channel)
        i = i + 1;
        InputWeights(Channel,i) = ANN(ANNindex).net.IW{1}(n);
    end
end

% Frequency magnitude and phase response
figure;
for n = 1:(DataLines)
    subplot(4,3,n);
    [h, f] = freqz(InputWeights(n,:),...
                   [1 zeros(1,120 - 1)],...
                   512,Parameters.SamplingFrequency);
    hAx = plotyy(f/1000,20*log10(abs(h)),f/1000,unwrap(angle(h)));
    ylabel(hAx(1),'Magnitude (dB)');
    ylabel(hAx(2),'Phase (rad)');
    %ylim(hAx(2),[-100 0]);
    ylim(hAx(1), [-60 60]);
    hAx(2).YTick = [-100 -75 -50 -25 0];
    xlim([0 0.5*Parameters.SamplingFrequency/1000]);
    xlabel('Frequency (kHz)');
    title(num2str(n));
end

% Phase delay
figure;
meanphi = zeros(size(DataLines));
for n = 1:DataLines
    [phi, w] = phasedelay(InputWeights(n,:),...
                         [1 zeros(1,120 - 1)],...
                         512);
    meanphi(n) = mean(phi(2:end));
    plot(w/pi,phi);
    hold on;
    labels{n} = ['Channel ' num2str(n)];
    xlim([0 1]);
    xlabel('Normalised Frequency (\times\pi rad/sample)');
    ylabel('Phase delay (samples)');
end
legend(labels);

% Plot average filter delay (phase delay) in each channel
figure;
plot(1:DataLines, meanphi, 'b*');
h_phase = lsline;   % Least-squares line
ChannelDelay = abs((h_phase.YData(1)-h_phase.YData(2))/...
               (h_phase.XData(1)-h_phase.XData(2)));
text(2, h_phase.YData(2)+3, ['Average inter-channel delay = ' num2str(ChannelDelay) ' samples']);
xlim([1, DataLines]);
xlabel('Channel');
ylabel('Average phase delay (samples)');

% % Group delay
% figure;
% for n = 1:(DataLines)
%     [gd, f] = grpdelay(InputWeights(n,:),...
%                          [1 zeros(1,120 - 1)],...
%                          512,Parameters.SamplingFrequency);
%     meangd(n) = mean(gd);
%     plot(f/1000,gd);
%     hold on;
%     labels{n} = ['Channel ' num2str(n)];
%     xlim([0 0.5*Parameters.SamplingFrequency/1000]);
%     xlabel('Frequency (kHz)');
%     ylabel('Group delay (samples)');
% end
% legend(labels);
% 
% % Plot average filter delay (group delay) in each channel
% figure;
% plot(1:DataLines, meangd, 'b*');
% h_group = lsline;
% ChannelDelay = abs((h_group.YData(1)-h_group.YData(2))/...
%                (h_group.XData(1)-h_group.XData(2)));
% text(2, h_group.YData(2)+3, ['Average inter-channel delay = ' num2str(ChannelDelay) ' samples']);
% xlim([1, DataLines]);
% xlabel('Channel');
% ylabel('Average group delay (samples)');

% Plot weights
figure;
MaxWeight = max(max(abs(InputWeights)));
for n = 1:(DataLines)
    subplot(DataLines,1,n)
    plot((0:120-1),...
         InputWeights(n,:));
    grid on;
    ylim([-MaxWeight MaxWeight]);
end
xlabel('Delayed samples (z^{-1})');

% % Plot FFT of weights
% figure;
% for n = 1:(DataLines)
%     subplot(DataLines,1,n)
%     Signal = InputWeights(n,:);
%     L = 2*floor(numel(Signal)/2);
%     Y = fft(Signal);
%     P2 = abs(Y/L);
%     P1 = P2(1:L/2+1);
%     P1(1:end) = 2*P1(1:end);
%     Frequency = Parameters.SamplingFrequency*(0:(L/2))/L/1000;
%     plot(Frequency,P1);
%     grid on;
%     MaxValue(n) = max(P1);
% end
% xlabel('Frequency (kHz)')
% for n = 1:(DataLines)
%     subplot(DataLines,1,n)
%     ylim([0 max(MaxValue)]);
% end

% % Plot phase response of filters
% figure;
% for n = 1:(DataLines)
%     subplot(DataLines,1,n)
%     Signal = InputWeights(n,:);
%     L = 2*floor(numel(Signal)/2);
%     Y = fft(Signal);
%     P = angle(Y);
%     phase = P(1:L/2+1);
%     Frequency = Parameters.SamplingFrequency*(0:(L/2))/L/1000;
%     plot(Frequency,phase);
%     grid on;
%     ylim([-pi pi]);
% end
% xlabel('Frequency (kHz)');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Investigate ANN outputs at each FIR filter system
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ANNindex = 1;

testnet = ANN(ANNindex).net;                % copy neural network
testnet.outputConnect = [1 0];              % connect output to first layer
testnet.biasConnect = [0; 0];
%testnet.layers{1}.transferFcn = 'purelin';  % 'purelin' for pre-sum

% Set the initial velocity set
Velocities = 10:10:120;
NumVelocities = numel(Velocities);

% Set up the time sequence for the signals at each velocity
InitTime = 0.001;
EndTime = 0.004;
Time = -InitTime:1/Parameters.SamplingFrequency:EndTime;
TimeLength = numel(Time);
SequenceTime = (0:(TimeLength * NumVelocities) - 1) ...
               / Parameters.SamplingFrequency;

% Run a test at each of the selected velocities
for VelocityIndex = 1:NumVelocities     
    
    % Retrieve the selected velocity
    Velocity = Velocities(VelocityIndex);
    
    % Get a data set of the approriate velocity
    Data = GetData(Parameters, Velocity, Time);
    
    % Insert the data into the input sequence
    InsertStart = 1 + (VelocityIndex - 1) * TimeLength;
    InsertEnd = InsertStart + TimeLength - 1;
    TripolarInputSequence (:,InsertStart:InsertEnd) = Data;
end

figure;
for ChannelIndex = 1:DataLines
    % Apply the AGC
    TripolarNorm = AgcSim(TripolarInputSequence);
    
    % Set every other channel to zero
    for i = 1:DataLines
        if i == ChannelIndex
            continue;
        end
        TripolarNorm(i,:) = zeros(1,numel(SequenceTime));
    end
    
    % Convert data format to be input to the ANN from a R-by-TS matrix
    % to a 1-by-TS cell array of R-by-1 vectors
    TripolarCellArray = con2seq(TripolarNorm(1:DataLines,:));
    
    % Run the ANN
    AnnOutputCellArray = sim(testnet,TripolarCellArray);
    
    % Convert ANN output back to a conventional array
    AnnOutput(ChannelIndex,:) = cell2mat(AnnOutputCellArray);
    title('FIR filter time response');
    subplot(DataLines,1,ChannelIndex);
    plot(SequenceTime, AnnOutput(ChannelIndex,:));
end
xlabel('Time (s)')

MaxValue = max(max(abs(AnnOutput)));
for n = 1:(DataLines)
    subplot(DataLines,1,n)
    ylim([-MaxValue MaxValue]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Investigate ANN outputs of first layer after summation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ANNindex = 1;
ChannelIndex = 1:DataLines;

testnet = ANN(ANNindex).net;                % copy neural network
testnet.outputConnect = [1 0];              % connect output to first layer
testnet.biasConnect = [1; 1];
%testnet.layers{1}.transferFcn = 'purelin';  % 'purelin' for pre-'tansig'

[TripolarInputSequence, Sequence] = signal_generator(Parameters, GetData);

% Apply the AGC
TripolarNorm = AgcSim(TripolarInputSequence);

% Convert data format to be input to the ANN from a R-by-TS matrix
% to a 1-by-TS cell array of R-by-1 vectors
TripolarCellArray = con2seq(TripolarNorm(ChannelIndex,:));

% Run the ANN
AnnOutputCellArray = sim(testnet,TripolarCellArray);

% Convert ANN output back to a conventional array
AnnOutput = cell2mat(AnnOutputCellArray);

figure;
subplot(2,1,1);
plot(SequenceTime, TripolarNorm(DataLines,:));
subplot(2,1,2);
plot(SequenceTime, AnnOutput);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Align filter weights
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ANNindex = 1:3
    MatchedVelocity = AnnVelocities(ANNindex);
    
    % Get weights
    InputWeights = zeros(DataLines,120);
    for Channel = 1:(DataLines)
        i = 0;
        for n = Channel:...
                (DataLines):...
                (numel(ANN(ANNindex).net.IW{1}) - (DataLines) + Channel)
            i = i + 1;
            InputWeights(Channel,i) = ANN(ANNindex).net.IW{1}(n);
        end
    end
   
    % Align weights
    figure;
    MaxWeight = max(max(abs(InputWeights)));
    DelaySumWeights = zeros(size(InputWeights(5,:)));
    InterElectrodeDelay = Parameters.ElectrodeSpacing*...
                          Parameters.SamplingFrequency/...
                          MatchedVelocity;
    for n = 1:DataLines
        Delay = round((n-5)*InterElectrodeDelay);
        DelayedWeights = circshift(InputWeights(n,1:end), Delay, 2);
        DelaySumWeights = DelaySumWeights + DelayedWeights;
        subplot(DataLines,1,n);
        plot((0:120-1), DelayedWeights);
        ylim([-MaxWeight MaxWeight]);
        grid on;
    end
    xlabel('Delayed samples (z^{-1})');
    
%     % Plot summed weights
%     figure;
%     MaxWeight = max(max(abs(DelaySumWeights)));
%     plot((0:120-1), DelaySumWeights);
%     grid on;
%     ylim([-MaxWeight MaxWeight]);
%     xlabel('Delayed samples (z^{-1})');

%     % Frequency magnitude and phase response
%     figure;
%     [h, f] = freqz(DelaySumOutput,...
%         [1 zeros(1,120 - 1)],...
%         500,Parameters.SamplingFrequency);
%     hAx = plotyy(f/1000,20*log10(abs(h)),f/1000,unwrap(angle(h)));
%     ylabel(hAx(1),'Magnitude (dB)');
%     ylabel(hAx(2),'Phase (rad)');
%     %ylim(hAx(2),[-100 0]);
%     %hAx(2).YTick = [-100 -75 -50 -25 0];
%     xlim([0 0.5*Parameters.SamplingFrequency/1000]);
%     xlabel('Frequency (kHz)');
%     title([num2str(MatchedVelocity) ' m/s']);
%     grid on;
    
%     testnet = ANN(ANNindex).net;                % copy neural network
%     testnet.outputConnect = [1 0];              % connect output to first layer
%     testnet.biasConnect = [1; 1];
%     testnet.layers{1}.transferFcn = 'purelin';  % 'purelin' for pre-'tansig'
%     testnet.input.size = 1;
%     testnet.IW{1} = DelaySumWeights;
%     
%     % Set the initial velocity set
%     Velocities = Parameters.StartTestVelocity:...
%                  Parameters.StepTestVelocity:...
%                  Parameters.EndTestVelocity;
%     NumVelocities = numel(Velocities);
%     
%     % Set up the time sequence for the signals at each velocity
%     InitTime = 0.001;
%     EndTime = 0.004;
%     Time = -InitTime:1/Parameters.SamplingFrequency:EndTime;
%     TimeLength = numel(Time);
%     SequenceTime = (0:(TimeLength * NumVelocities) - 1) ...
%         / Parameters.SamplingFrequency;
%     
%     % Run a test at each of the selected velocities
%     for VelocityIndex = 1:NumVelocities
%         
%         % Retrieve the selected velocity
%         Velocity = Velocities(VelocityIndex);
%         
%         % Get a data set of the approriate velocity
%         TripolarData = GetTriPolar(Parameters, Velocity, Time);
%         
%         % Insert the data into the input sequence
%         InsertStart = 1 + (VelocityIndex - 1) * TimeLength;
%         InsertEnd = InsertStart + TimeLength - 1;
%         TripolarInputSequence (:,InsertStart:InsertEnd) = TripolarData;
%     end
%     
%     % Apply the AGC
%     TripolarNorm = AgcSim(TripolarInputSequence);
%     
%     % Delay and sum
%     DelaySumOutput = TripolarNorm(DataLines,:);
%     InterElectrodeDelay = Parameters.ElectrodeSpacing/MatchedVelocity;
%     for n = 1:Parameters.Electrodes-3
%         Delay = round((DataLines-n)*...
%                       InterElectrodeDelay*...
%                       Parameters.SamplingFrequency);
%         DelayedSignal = circshift(TripolarNorm(n,1:end), Delay, 2);
%         DelaySumOutput = DelaySumOutput + DelayedSignal;
%     end
%     
%     % Convert data format to be input to the ANN from a R-by-TS matrix
%     % to a 1-by-TS cell array of R-by-1 vectors
%     TripolarCellArray = con2seq(DelaySumOutput);
%     
%     % Run the ANN
%     AnnOutputCellArray = sim(testnet,TripolarCellArray);
%     
%     % Convert ANN output back to a conventional array
%     AnnOutput = cell2mat(AnnOutputCellArray);
%     
%     figure;
%     subplot(2,1,1);
%     plot(SequenceTime, DelaySumOutput);
%     subplot(2,1,2)
%     plot(SequenceTime, AnnOutput);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SNR estimation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Parameters = struct (                 ...
   'Electrodes',              11    , ...
   'ElectrodeSpacing',      3e-3    , ...
   'SamplingFrequency',      1e5    , ...
   'ActionPotentialType',    'tmap2', ...
   'StartTestVelocity',       10    , ...
   'StepTestVelocity',         5    , ...
   'EndTestVelocity',        120    , ...
   'MatchRepeats',            10    , ...
   'NoiseLevel',           1e-21    );

% Create a noisy signal with a range of velocities
[TriPolarSignal, ~] = signal_generator(Parameters, GetData);

% Create a noiseless signal with a range of velocities
Parameters.NoiseLevel = 0;
[NoiselessTriPolarSignal, SequenceTime] = signal_generator(Parameters, GetData);

% Get noise record
Noise = TriPolarSignal - NoiselessTriPolarSignal;

% SNR in each channel
for i = 1:DataLines
    SignalNoiseRatio(i) = snr(NoiselessTriPolarSignal(i,:), Noise(i,:));
end

% Average SNR across channels
SNR = mag2db(mean(db2mag(SignalNoiseRatio)));

% Plot signal and noise
figure;
hold on;
plot(SequenceTime,NoiselessTriPolarSignal(1,:));
plot(SequenceTime,Noise(1,:));
text(0.1*max(SequenceTime), max(abs(NoiselessTriPolarSignal(1,:))),...
     ['SNR = ' num2str(SNR) ' dB']);
xlim([0 0.005*NumVelocities]);
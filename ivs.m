function ivs()
% Plot the intrinsic velocity spectrum for comparison with the ANNs
%
% A Vlissidis

%% Delay and sum method
Parameters = struct (                 ...
   'Electrodes',              11    , ...
   'ElectrodeSpacing',      3e-3    , ...
   'SamplingFrequency',      1e5    , ...
   'ActionPotentialType',    'tmap2', ...
   'StartTestVelocity',       10    , ...
   'StepTestVelocity',       0.1    , ...
   'EndTestVelocity',        160    , ...
   'APType',              'TriPolar', ...
   'NoiseLevel',               0    );

if Parameters.APType == 'UniPolar'
    GetData = @GetUniPolar;
    DataLines = Parameters.Electrodes;
elseif Parameters.APType == 'TriPolar'
    GetData = @GetTriPolar;
    DataLines = Parameters.Electrodes - 2;
end

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

end


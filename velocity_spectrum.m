function velocity_spectrum(ANN, Parameters)
% Plot the velocity spectrum of the trainned model
%
% A Vlissidis

%MatchedVelocities = [80, 90, 120];
MatchedVelocities = [10, 20, 30];

% Velocities to test at
Velocities = Parameters.StartTestVelocity:...
             Parameters.StepTestVelocity:...
             Parameters.EndTestVelocity;
NumVelocities = numel(Velocities);

% Set up the time sequence for the signals
InitTime = 0.001;
EndTime = 0.004;
Time = -InitTime:1/Parameters.SamplingFrequency:EndTime;

NumAnns = numel(ANN);
PeakBuiltIn = zeros(NumVelocities,NumAnns);
PeakCustom  = zeros(NumVelocities,NumAnns);

% Handle to set action potential type
if strcmp(Parameters.APType, 'UniPolar')
    GetData = @GetUniPolar;
    DataLines = Parameters.Electrodes;
elseif strcmp(Parameters.APType, 'TriPolar')
    GetData = @GetTriPolar;
    DataLines = Parameters.Electrodes - 2;
elseif strcmp(Parameters.APType, 'BiPolar')
    GetData = @GetBiPolar;
    DataLines = floor(Parameters.Electrodes/2);
end

% Set up an array for the input tripoles so that we can do an overall AGC
Data = zeros(NumVelocities , DataLines , numel(Time));

% Run a test at each of the selected velocities
for VelocityIndex = 1:NumVelocities     
    
    % Retrieve the selected velocity
    Velocity = Velocities(VelocityIndex);
    
    % Get a data set of the approriate velocity using default electrode
    % parameters and no noise for a theoretical tmap.
    Data(VelocityIndex,:,:) = GetData(Parameters, Velocity, Time);
end

% Apply the AGC
TripolarNorm = AgcSim(Data);

for VelocityIndex = 1:NumVelocities          
    % Convert data format to be input to the ANN from a R-by-TS matrix
    % to a 1-by-TS cell array of R-by-1 vectors
    TripolarCellArray = con2seq(reshape(TripolarNorm(VelocityIndex,:,:),...
        size(TripolarNorm,2),size(TripolarNorm,3))); 

    % Apply data of the selected velocity to each ANN
    for AnnIndex=1:NumAnns
        % Run the ANN
        AnnOutputCellArray = sim(ANN(AnnIndex).net,TripolarCellArray);
         
        % Convert ANN output back to a conventional array
        AnnOutput = cell2mat(AnnOutputCellArray);
        
        % Get the peak reponse
        PeakBuiltIn(VelocityIndex,AnnIndex) = max(AnnOutput);
        PeakCustom(VelocityIndex,AnnIndex) = max(abs(AnnOutput));
    end 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Use the peak responses to plot an IVS for each ANN 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot the peaks
figure;
title('Peak responses using built in ANN tools');
xlabel('Velocity (m/s)');
ylabel('Magnitude of response');
hold all;

% plot each of the selected ANNs
for AnnIndex=1:NumAnns
    
    % Calculate velocity selectivity
    v3 = Velocities(PeakBuiltIn(:,AnnIndex) >= (...
              PeakBuiltIn(Velocities == MatchedVelocities(AnnIndex),AnnIndex)...
              /sqrt(2)));
    BW = max(v3) - min(v3);
    
    MatchedVelocities(AnnIndex);
    Qv = MatchedVelocities(AnnIndex) / BW;
    
    % Plot the peak reponse
    plot (Velocities,PeakBuiltIn(:,AnnIndex),...
          'DisplayName',[num2str(MatchedVelocities(AnnIndex)) ' m/s']);
    
    text(MatchedVelocities(AnnIndex),...
         max(PeakBuiltIn(:,AnnIndex)), ['Q_v = ' num2str(Qv)]);
end

% Show the legend
legend('show');

end


function ANN = DesignAnn (AnnSelections, Parameters)
% Design ANNs using the built in training systems to detect APs at
% particular velocities.
%
% C.T.Clarke based on the work of Assad al Shueli
% Edited by L F Tiong 06/05/2016
%
% Use a parameter structure like this:
%
% Parameters = struct (                ...
%   'Electrodes',              11    , ...
%   'ElectrodeSpacing',         0.003, ...
%   'SamplingFrequency',   100000    , ...
%   'ActionPotentialType', 'long'    , ...
%   'StartTestVelocity',        1    , ...
%   'StepTestVelocity',         1    , ...
%   'EndTestVelocity',        100    , ...
%   'NoiseLevel',               0.01 );
%
% Other parameters may be included but will be ignored

% Load all ANNs which have been stored as array called 
% ANN(i).net where i is the ANN index that includes all 
% the values in AnnSelections below.

% Get the number of ANNs to train
NumAnns = numel(AnnSelections);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup required variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set up the time sequence for the signals at each velocity
InitTime = 0.001;
EndTime = 0.004;
Time = -InitTime:1/Parameters.SamplingFrequency:EndTime;
TimeLength = numel(Time);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Design each ANN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for AnnIndex=1:NumAnns

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Create input signals for training
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Retrieve the matched velocity of the selected ANN
    MatchedVelocity = AnnSelections(AnnIndex);

    % Set the initial Velocity set including the repeated matched set at
    % the end that match the ANN's target velocity
    Velocities = cat(2,Parameters.StartTestVelocity: ...
        Parameters.StepTestVelocity:Parameters.EndTestVelocity, ...
        MatchedVelocity*ones(1,Parameters.MatchRepeats));
    NumVelocities = numel(Velocities);
    
    % Set up the array to hold the input training data
    TrainingInputSequence = zeros(Parameters.Electrodes - 2, ...
                              NumVelocities * TimeLength);
                          
    % Create an input sequence for the ANN
    for VelocityIndex = 1:NumVelocities     

        % Retrieve the selected velocity
        Velocity = Velocities(VelocityIndex);

        % Get a data set of the appropriate velocity
        TripolarData = GetTriPolar(Parameters, Velocity, Time);

        % Insert the data into the input sequence
        InsertStart = 1 + (VelocityIndex - 1) * TimeLength;
        InsertEnd = InsertStart + TimeLength - 1;
        TrainingInputSequence (:,InsertStart:InsertEnd) = TripolarData;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Create an output target signal for training
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Clear the target
    TrainingTargetSequence = zeros(1,NumVelocities*TimeLength);

    % Add a pulse for each velocity that is the same as the matched
    % velocity 
    VelocityIndices = 1:NumVelocities;
    for VelocityIndex = VelocityIndices(Velocities == MatchedVelocity)    
             
        % Get the insertion point for the data in the input sequence
        InsertStart = 1 + (VelocityIndex - 1) * TimeLength;
        InsertEnd = InsertStart + TimeLength - 1;

        % Get a data set of the appropriate velocity with no noise
        NoiseLessParameters = Parameters;
        NoiseLessParameters.NoiseLevel = 0;
        TripolarData = GetTriPolar(NoiseLessParameters, MatchedVelocity, Time);

        % Get the final channel
        FinalChannel = TripolarData(Parameters.Electrodes - 2,:);

        % Normalise it 
        FinalChannel = AgcSim(FinalChannel);
        
        switch lower(Parameters.TargetOutput)
            case 'pulse'
                % Create a pulse at the 3dB points
                TargetPulse = zeros(size(FinalChannel));
                TargetPulse(FinalChannel >= (max(FinalChannel)/sqrt(2))) = 1;
                
                % Apply it to the target sequence
                TrainingTargetSequence(InsertStart:InsertEnd) = TargetPulse;
            case 'ap'
                % Target sequence is action potential with matched velocity
                TrainingTargetSequence(InsertStart:InsertEnd) = FinalChannel;
            otherwise
                disp('Unknown output training target type.')
        end;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Create input signals for validation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    valParameters = Parameters;
    valParameters.StepTestVelocity = 5;
    valParameters.NoiseLevel = 2e-21;
    valParameters.MatchRepeats = 5;
    
    % Set the validation Velocity set including the repeated matched set at
    % the end that match the ANN's target velocity
    Velocities = cat(2,valParameters.StartTestVelocity: ...
        valParameters.StepTestVelocity:valParameters.EndTestVelocity, ...
        MatchedVelocity*ones(1,valParameters.MatchRepeats));
    NumVelocities = numel(Velocities);
    
    % Set up the array to hold the input training data
    ValidationInputSequence = zeros(valParameters.Electrodes - 2, ...
                              NumVelocities * TimeLength);
                          
    % Create an input sequence for the ANN
    for VelocityIndex = 1:NumVelocities     

        % Retrieve the selected velocity
        Velocity = Velocities(VelocityIndex);

        % Get a data set of the appropriate velocity
        TripolarData = GetTriPolar(valParameters, Velocity, Time);

        % Insert the data into the input sequence
        InsertStart = 1 + (VelocityIndex - 1) * TimeLength;
        InsertEnd = InsertStart + TimeLength - 1;
        ValidationInputSequence (:,InsertStart:InsertEnd) = TripolarData;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Create an output target signal for validation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Clear the target
    ValidationTargetSequence = zeros(1,NumVelocities*TimeLength);

    % Add a pulse for each velocity that is the same as the matched
    % velocity 
    VelocityIndices = 1:NumVelocities;
    for VelocityIndex = VelocityIndices(Velocities == MatchedVelocity)    
             
        % Get the insertion point for the data in the input sequence
        InsertStart = 1 + (VelocityIndex - 1) * TimeLength;
        InsertEnd = InsertStart + TimeLength - 1;

        % Get a data set of the appropriate velocity with no noise
        valParameters.NoiseLevel = 0;
        TripolarData = GetTriPolar(valParameters, MatchedVelocity, Time);

        % Get the final channel
        FinalChannel = TripolarData(valParameters.Electrodes - 2,:);

        % Normalise it 
        FinalChannel = AgcSim(FinalChannel);
        
        switch lower(valParameters.TargetOutput)
            case 'pulse'
                % Create a pulse at the 3dB points
                TargetPulse = zeros(size(FinalChannel));
                TargetPulse(FinalChannel >= (max(FinalChannel)/sqrt(2))) = 1;
                
                % Apply it to the target sequence
                ValidationTargetSequence(InsertStart:InsertEnd) = TargetPulse;
            case 'ap'
                % Target sequence is action potential with matched velocity
                ValidationTargetSequence(InsertStart:InsertEnd) = FinalChannel;
            otherwise
                disp('Unknown output training target type.')
        end;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Run the ANN training
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Normalise the data into the range 1 to -1 for each channel separately
    TripolarNorm = AgcSim([TrainingInputSequence, ValidationInputSequence]);
    %TripolarNorm = AgcSim(TrainingInputSequence);
    
    TargetSequence = [TrainingTargetSequence, ValidationTargetSequence];
    %TargetSequence = TrainingTargetSequence; 
    
    % Assad's ANN training setup
    p = con2seq(TripolarNorm);      % Cell array of input vectors
    t = con2seq(TargetSequence);    % Cell array of target vectors
    
    for i = 1:numel(Parameters.AnnLength)
        d{i} = 0:Parameters.AnnLength(i) - 1;    % Delay vector for ith layer
    end
    
    %net = newdtdnn(p,t,1,d);         % Obsolete function
    net = distdelaynet(d, Parameters.HiddenLayerSize, Parameters.TrainingMethod);
    
    net.layers{1:end-1}.transferFcn = Parameters.ActivationFcn1;
    net.layers{end}.transferFcn = Parameters.ActivationFcn2;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Levenberg-Marquardt backpropagation default parameters
    % net.trainParam.epochs = 1000;         % Maximum number of epochs to train
    % net.trainParam.goal = 0;              % Performance goal
    % net.trainParam.max_fail = 6;          % Maximum validation failures
    % net.trainParam.min_grad = 1e-7;       % Minimum performance gradient
    % net.trainParam.mu = 0.001;            % Initial mu
    % net.trainParam.mu_dec = 0.1;          % mu decrease factor
    % net.trainParam.mu_inc = 10;           % mu increase factor
    % net.trainParam.mu_max = 1e10;         % Maximum mu
    % net.trainParam.show = 25;             % Epochs between displays (NaN for no displays)
    % net.trainParam.showCommandLine = 0;   % Generate command-line output
    % net.trainParam.showWindow = 1;        % Show training GUI
    % net.trainParam.time = Inf;            % Maximum time to train in seconds
    % net.trainParam.mem_reduc = 1;         % Reduce memory and speed to calculate the Jacobian jX
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    net.trainParam.max_fail = 5;
    net.divideFcn = 'divideind'; 
    net.divideParam.trainInd = 1:size(TrainingInputSequence,2);
    net.divideParam.valInd = (1:size(ValidationTargetSequence,2)) + ...
                             size(TrainingInputSequence,2);
    
    % These are Assad's training criteria so they have been kept
    net.trainParam.epochs = Parameters.MaxIterations;
    net.trainParam.goal = Parameters.MinError;
    %net.divideFcn = 'dividetrain';        % Allocate all data for training
    %net.trainParam.lr = 0.00005;          % trainlm does not have a learning rate?
    
    % Do the ANN training
    net = train(net,p,t);
    
    % Store the ANN in an array 
    ANN(AnnIndex).net = net;
end
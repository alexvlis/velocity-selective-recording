function time_analysis (channels, ANNindex, opt, rawdata)
% Analyse the time response of the neural network
%
% A Vlissidis

    ANN = evalin('base', 'ANN');
    Parameters = evalin('base', 'Parameters');
    
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

    %% Create the input data
    % Set the initial velocity set
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

    for VelocityIndex = 1:NumVelocities     

        % Retrieve the selected velocity
        Velocity = Velocities(VelocityIndex);

        % Get a data set of the approriate velocity
        Data = GetData(Parameters, Velocity, Time);

        % Insert the data into the input sequence
        InsertStart = 1 + (VelocityIndex - 1) * TimeLength;
        InsertEnd = InsertStart + TimeLength - 1;
        InputSequence (:,InsertStart:InsertEnd) = Data;
    end
    
    if strcmp(opt, 'real')
       InputSequence = rawdata(1:7000, 2:end)';
       SequenceTime = rawdata(1:7000, 1)';
    end
    
    %% Plot each FIR filter output individually
    testnet = ANN(ANNindex).net;                % copy neural network
    testnet.outputConnect = [1 0];              % connect output to first layer
    testnet.biasConnect = [0; 0];
    testnet.layers{1}.transferFcn = 'purelin';
    
    figure;
    hold on
    % Pre-allocate for speed
    AnnOutput = zeros(DataLines, numel(SequenceTime));
    for ChannelIndex = 1:DataLines
        % Apply the AGC
        NormInput = AgcSim(InputSequence);

        % Set every other channel to zero
        NormInput(1:end ~= ChannelIndex, :) = zeros(DataLines-1, numel(SequenceTime));
        
        % Convert data format to be input to the ANN from a R-by-TS matrix
        % to a 1-by-TS cell array of R-by-1 vectors
        InputCellArray = con2seq(NormInput);

        % Run the ANN
        AnnOutputCellArray = sim(testnet, InputCellArray);
        
        % Convert ANN output back to a conventional array
        AnnOutput(ChannelIndex,:) = cell2mat(AnnOutputCellArray);
        
        % Plot channel time response
        %subplot(DataLines, 1, ChannelIndex);
        plot(SequenceTime, AnnOutput(ChannelIndex, :));
        if ChannelIndex == 1
            title('Individual FIR filter outputs');
        end
    end
    xlabel('Time (ms)');
    hold off
    %legend('1', '2', '3')
    %% Plot the summed output of the specified channels
    figure;
    % Sum the specified channels
    out = sum(AnnOutput(channels, :), 1);
    
    NormInput = AgcSim(InputSequence);
    InputCellArray = con2seq(NormInput);
    
    AnnOutputCellArray = sim(testnet, InputCellArray);
    AnnOutput = cell2mat(AnnOutputCellArray);
    sum(AnnOutput - out)
    
    plot(SequenceTime, out);
    xlabel('Time (ms)');
    title(strcat('Summed output from FIR Channels'));
    
    %% Plot the input
    if strcmp(opt, 'real')
        figure;
        for i = 2:6
            subplot(5, 1, i-1);
            plot(rawdata(1:7000, 1),rawdata(1:7000, i));
        end
        xlabel('Time (s)');
    else    
        figure;
        input = cell2mat(InputCellArray);
        plot(SequenceTime, input(DataLines, :));
        xlabel('Time (ms)');
        title('Input Sequence');
    end
    
    %% Plot summed output with bias
    testnet = ANN(ANNindex).net;
    testnet.outputConnect = [1 0];
    testnet.layers{1}.TransferFcn = 'purelin';
    
    Norm = AgcSim(InputSequence);
    InputCellArray = con2seq(Norm);
    AnnOutputCellArray = sim(testnet, InputCellArray);
    AnnOutput = cell2mat(AnnOutputCellArray);
    figure;
    plot(SequenceTime, AnnOutput);
    title('Summed output with bias');
    xlabel('Time (ms)');
    
    %% Plot summed output with bias and sigmoid function
    testnet = ANN(ANNindex).net;
    testnet.outputConnect = [1 0];
    testnet.layers{1}.transferFcn = 'tansig';
    
    AnnOutputCellArray = sim(testnet, InputCellArray);
    AnnOutput = cell2mat(AnnOutputCellArray);
    figure;
    plot(SequenceTime, AnnOutput);
    title('Summed output with bias and sigmoid');
    xlabel('Time (ms)');
   
    %% Revert back to normal and plot global output
    testnet = ANN(ANNindex).net;
    testnet.outputConnect = [0 1];
    testnet.layers{1}.transferFcn = 'tansig';
    
    AnnOutputCellArray = sim(testnet, InputCellArray);
    AnnOutput = cell2mat(AnnOutputCellArray);
    figure;
    plot(SequenceTime, AnnOutput);
    title('Global Output');
    xlabel('Time (ms)');
end
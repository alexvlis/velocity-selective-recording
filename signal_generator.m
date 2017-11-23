function [ signal, SequenceTime ] = signal_generator( Parameters, GetData )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

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
    signal(:,InsertStart:InsertEnd) = Data;
end

end


function TriPolarSignal = GetTriPolar(Parameters, ActionPotentialVelocity, Time)
% Create a Tripolar action potential
%
% C T Clarke
% Edited by L F Tiong 06/05/2016
%
% Use a parameter structure like this:
% Parameters = struct (                ...
%   'Electrodes',              11    , ...
%   'ElectrodeSpacing',         0.003, ...
%   'SamplingFrequency',   100000    , ...
%   'ActionPotentialType', 'long'    , ...
%   'StartTestVelocity',        1    , ...
%   'StepTestVelocity',         1    , ...
%   'EndTestVelocity',        100    , ...
%   'NoiseLevel',               0.01 );
% ElectrodeSpacing is not used here but is needed by GetUniPolar
% Other parameters may be included but will be ignored

% Work out how many tripoles there will be
Tripoles = Parameters.Electrodes - 2;

% Work out how long the array should be
SequenceLength = max(size(Time));

% Standard Unipolar signal
UniPolarSignal = GetUniPolar(Parameters, ActionPotentialVelocity, Time);

% Resultant Tripole
TriPolarSignal = zeros(Tripoles , SequenceLength);
for Count = 1 : Tripoles
    TriPolarSignal(Count,:) = UniPolarSignal(Count,:) - 2 * UniPolarSignal(Count + 1,:) + UniPolarSignal(Count + 2,:); 
end


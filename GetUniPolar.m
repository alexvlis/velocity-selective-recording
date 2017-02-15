function UniPolarSignal = GetUniPolar(Parameters, ActionPotentialVelocity, Time)
% Create a Unipolar action potential
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
% Other parameters may be included but will be ignored

% Establish implied cuff parameters
if (Parameters.Electrodes > 1)
    ElectrodeDelay = Parameters.ElectrodeSpacing / ActionPotentialVelocity; 
else
    ElectrodeDelay = 0;
end

% Get the length of sequence to produce
SequenceLength = max(size(Time));


A = [2.2e7,0.47e9,2.6e1,4.08e-3,7.44e-11];
B = [3.6e3,  1e4, 1.5e4,  1.5e4,    1e4];
n = [3,        3,     1,      1,      3];
%Use 1 for long, 2 for short, 3 for new
                   
% Use on of the standard formulae for the action potential    
switch lower(Parameters.ActionPotentialType)
    case 'long'
        TMAP = 1;
    case 'short'
        TMAP = 2;
    case 'new'
        TMAP = 3;
    case 'tmap1'
        TMAP = 4;
    case 'tmap2'
        TMAP = 5;
    otherwise
        disp('Unknown Action potential type.')
end;

% Standard Unipolar signal (ActionPotentialVelocity^2 * )
UniPolarSignal = zeros(Parameters.Electrodes,SequenceLength);
for Count = 1 : Parameters.Electrodes
    UniPolarSignal(Count,:) = ActionPotentialVelocity^2 * max(0,A(TMAP)*((Time - (Count-1) * ElectrodeDelay).^n(TMAP)).*(exp((-B(TMAP))*(Time - (Count-1) * ElectrodeDelay))));
end;

% Work out the approximate RMS value of the signal
% max/root(2) is used as this negates the effect of an overly long sequence
% with a short AP.
%SignalRMS = max(max(UniPolarSignal)) / sqrt(2);

% Constant noise level for all AP amplitudes (SignalRMS * )

% Add normally distributed noise
UniPolarSignal = UniPolarSignal + Parameters.NoiseLevel * randn(Parameters.Electrodes,SequenceLength);


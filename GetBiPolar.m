function BiPolarSignal = GetBiPolar(Parameters, ActionPotentialVelocity, Time)
% Create a BiPolar action potential
%
% A Vlissidis

BiPoles = floor(Parameters.Electrodes/2);
SequenceLength = max(size(Time));

UniPolarSignal = GetUniPolar(Parameters, ActionPotentialVelocity, Time);

idx = 1;
BiPolarSignal = zeros(BiPoles, SequenceLength);
for i = 1:2:Parameters.Electrodes - 1
    BiPolarSignal(idx, :) = UniPolarSignal(i, :) - UniPolarSignal(i+1, :);
    idx = idx + 1;
end

end
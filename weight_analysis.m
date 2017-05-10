function weight_analysis(ann, Parameters)
% Analyse the weight patterns
%  
% A. Vlissidis

weights = zeros(Parameters.channels, Parameters.AnnLength(1));
% Plot the weights of the FIR filters
figure;
for Channel = 1:Parameters.channels
    i = 0;
    for n = Channel:Parameters.channels:...
            (numel(ann.IW{1}) - (Parameters.channels) + Channel)
        i = i + 1;
        weights(Channel,i) = ann.IW{1}(n);
    end
    subplot(Parameters.channels, 1, Channel);
    plot(weights(Channel, :));
end

% Plot the phase response of the FIR filters
figure;
for Channel = 1:Parameters.channels
    subplot(Parameters.channels, 1, Channel);
    phasez(weights(Channel, :), 1, Parameters.AnnLength(1));
end

figure;
[phi, f] = phasedelay(weights(1, :), 1, 512, 1e5);
plot(f/1000, phi);

mean(weights(1, :))
mean(weights(2, :))
mean(weights(3, :))
% mean(weights(4, :))
% mean(weights(5, :))
end
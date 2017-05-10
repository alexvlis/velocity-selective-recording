function h = sim_fir(b, c, Parameters, x)
%SIM_FIR Summary of this function goes here
%   Detailed explanation goes here

a = 1;
y = zeros(size(x)); % Allocate the space for the output

for channel = 1:Parameters.channels
    % Apply centroid
    x(channel, :) = filter(b(channel, :), a, x(channel, :));
    
    if channel < Parameters.channels
        % Apply delay
        y(channel, :) = filter(c(channel, :), a, x(channel, :));
    else
        y(channel, :) = 1.813 * x(channel, :); % scale last channel
    end

end

% Sum the output of each channel
h = sum(y);
end


function Result = AgcSim(Data)
% Approximate the effect of individual channel normalisation and DC removal
% on data. The output data has mean value and either a maximum value of 1
% or a minimum value of -1. Can be used on arrays of multi-channel data. 
% It is assumed that the final dimension is time and the dimension before
% that is the channels. If there is a third dimension, then dimension 1 is
% the seqeunce set. Normalisation is applied to the the all the data for a
% channel creating the efect of feeding the array through the channels as a
% continuous stream.
% C. T. Clarke
% Edited by L F Tiong 06/05/2016

DataNoDc = zeros (size(Data));
Result = zeros (size(Data));

% Normalise the data into the range 1 to -1 for each channel separately
for Channel = 1:size(Data,ndims(Data)-1)
    
    % DC removal
    if (ndims(Data) == 2) 
       DataNoDc(Channel,:) = Data(Channel,:) - mean(Data(Channel,:));
    else
       DataNoDc(:,Channel,:) = Data(:,Channel,:) - mean(mean(Data(:,Channel,:)));
    end
    
    % Normalize range to +/-1.0
    if (ndims(Data) == 2) 
        Result(Channel,:) = DataNoDc(Channel,:) / max(abs(DataNoDc(Channel,:)));
    else
        Result(:,Channel,:) = DataNoDc(:,Channel,:) / max(max(abs(DataNoDc(:,Channel,:))));
    end
end

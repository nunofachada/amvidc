function data_norm = removeMean(data)
% REMOVEMEAN Removes matrix mean vector from each row.
%
% data_norm = REMOVEMEAN(data)
%
% Parameters:
%       data - data from which to remove mean from each row
% Output:
%  data_norm - data with mean removed from each row
%

%  N. Fachada
%  Instituto Superior Técnico, Lisboa, Portugal

% Get mean
meanVector = mean(data);

% Pre-allocate
data_norm = zeros(size(data));

% Remove mean, row by row
for i=1:size(data, 1)
    data_norm(i, :) = data(i, :) - meanVector; 
end;

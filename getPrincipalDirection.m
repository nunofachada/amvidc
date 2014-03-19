function direction = getPrincipalDirection(data, idx, numClusts, method)
% GETPRINCIPALDIRECTION Return the principal direction of the given data, 
% possibly differentiating and weighting different clusters.
%
% direction = GETPRINCIPALDIRECTION(data, idx, numClusts, method)
%
% Parameters:
%        data - m x n, m samples, n dimensions
%         idx - m x 1, maps each sample with a cluster
%   numClusts - number of different clusters
%      method - 'pca' or 'svd' ('pca' obtains eigenvalues/weights using
%               the data covariance matrix, while 'svd' obtains them
%               directly from the data)
% Output:
%   direction - n x 1, normalized cluster-weighted principal direction
%

%  N. Fachada
%  Instituto Superior Técnico, Lisboa, Portugal

% Pre-alocate space for global direction data
globalDirectionVectors = zeros(size(data, 2), numClusts);
globalDirectionWeights = zeros(1, numClusts);

% Determine principal directions and weights of all clusters
for i=1:numClusts
    
    % Get data in cluster i
    currData = data(find(idx == i), :);
    
    % Cluster must have more than one sample to yield a direction
    if size(currData, 1) > 1
        
        % Get principal direction of current data using the specified method
        if strcmp(method, 'svd') 

            % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
            % Use SVD: eigenvalues are taken directly from the data matrix  %
            % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

            % SVD requires that we specifically remove the mean
            currData = removeMean(currData);

            % Doing this would yield same results as PCA method
            %currData = cov(currData);

            % Perform SVD
            [~, Ss, Vs] = svds(currData, 1);

        elseif strcmp(method, 'pca')

            % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
            % Use SVD: eigenvalues are taken directly from the data matrix  %
            % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

            % Perform PCA
            [Vs, ~, Ss] = princomp(currData, 'econ');
            Vs = Vs(:, 1);

        else
            error('Unknown method for obtaining principal direction: "%s"', method);
        end;
    else
        % Cluster only has one sample, give it zero weight
        Vs = zeros(size(data, 2), 1);
        Ss = 0;
    end;

    % Keep PD vector and respective weight
    globalDirectionVectors(:, i) = Vs;
    globalDirectionWeights(i) = abs(Ss(1));
    
end;

% Determine weighted global direction
direction = weightedDirection(globalDirectionVectors, globalDirectionWeights);

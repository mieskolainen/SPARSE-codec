% K-means quantization/clustering of the data vectors
% ------------------------------------------------------------------------
% Cluster into K classes using Euclidean L_2 distance
%
% Input:         X  : d x N matrix of data vectors
%                K  : number of clusters
%            iters  : maximum number of iterations (default = inf)
%
% Output:      cent : dxK vector of quantization points
%             clust : index of quantization point for each data vector
%               err : total squared error of quantized data
% 
% Mikael Mieskolainen, mikael.mieskolainen@tut.fi, 2011
% ------------------------------------------------------------------------

function [cent, clust, err] = K_means(X, K, iters)

if ( nargin < 3),
    iters = inf;
end

[D N] = size(X);

if (N < K)
    error(['K-means():: Error, we have a trivial solution, N < K ( ' ...
           int2str(N) ' < ' int2str(K) ')']);
end

iter = 1;
done = 0;

% Initialize centroids
[~,ind] = sort(rand(1,N));
cent = X( :,ind(1:K));

dist = zeros( K, N);
X2 = sum(sum(X.^2));

while ((~done ) && ( iter < iters))
    
    % Find distances to centroids
    for i = 1:K 
        if ( D > 1)
            dist( i,:) = sum((X - cent(:,i)*ones(1,N)).^2);
        else
            dist( i,:) = (X - cent(:,i)*ones(1,N)).^2;
        end
    end
    
    % Classify
    [~, clust] = min(dist);
    mindist = dist((0:N - 1)*K + clust); % find smallest distances
    
    % check if done
    if (iter == 1)
        prevclust = zeros(size(clust));
    end

    if (clust == prevclust) % We are done
        done = 1;
    else
        numchanges = sum(clust ~= prevclust);
        prevclust = clust;
    end
    
    err = sum( mindist);
    if (rem(iter, 10) == 1)
        fprintf('K-Means():: Iteration %d  SNR: %0.3f, Changes %d\n', ...
                iter, 10*log10( X2/err), numchanges);
    end
    
    % Calculate new centroids as means of the clusters
    for i = 1:K
        if ( sum( clust == i) == 0) % Empty cluster, pick new one
            cent(:, i) = X(:,ceil( rand*N));
        else
            cent(:, i) = mean(X(:, clust == i).').';
        end
    end
    
    iter = iter + 1;
end
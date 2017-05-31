% K-SVD sparse overcomplete dictionary learning algorithm
% ------------------------------------------------------------------------
% 
% Input:      X  =  Signal matrix (d x M), where d = dimension and
%                   M is the number of signal vectors
%             K  =  Number of Dictionary elements (atoms) to be trained 
%                   where usually K << M
%         iters  =  Number of iterations, for example 40
%          init  =  'Self', 'Gaussian'
%             L  =  Maximum number of coefficients in pursuit algorithm
%
% Output:     D  =  Trained dictionary (d x K)
%             A  =  Coefficient matrix (K x M)
%
% ------------------------------------------------------------------------
% M. Aharon, M. Elad, and A.M. Bruckstein. 
% The K-SVD: An Algorithm for Designing of Overcomplete Dictionaries for
% Sparse Representation. IEEE Trans. On Signal Processing, Vol. 54,
% no. 11, pp. 4311-4322, November 2006. 
%
% Mikael Mieskolainen, mikael.mieskolainen@tut.fi, 2011
% ------------------------------------------------------------------------

function [D, A] = ksvd(X, K, iters, init, L)

% ------------------------------------------------------------------------
% 0. Initialization part

[d,M] = size(X);

% If we have a trivial dictionary
if (M < K)
    fprintf('We have a trivial dictionary (M < K) (%d < %d)\n', M, K);
    D = normD(X);
    return;
end

% Random Initialization of dictionary
if (strcmp(init, 'Self'))
    randcols = randperm(M);
    D = X(:,randcols(1:K));
elseif (strcmp(init, 'Gaussian'))
    D = randn(d,K); 
end

A = sparse(K, M);  % Initialize the sparse coefficient matrix
D = normD(D);      % Normalize the dictionary

fprintf('\nStarting K-SVD algorithm... \n');

% Through every iteration
for i = 1:iters
    
    % --------------------------------------------------------------------
    % 1. SPARSE CODING
    
    tic;
    
    % Calculate the coefficient matrix using orthogonal matching pursuit
    A = OMP(D, X, L);
    
    % --------------------------------------------------------------------
    % 2. DICTIONARY UPDATE
    % Codebook update stage for j = 1,2, ... K with permutated order
    
    permutations = randperm(size(D,2));
    
    for j = permutations
        
        % Find a better atom vector d. Note: update A also at same time!
        [d, A] = getd(X, D, j, A);

        % Replace the previous atom vector with the new d
        D(:,j) = d;
    end
    
    time = toc;
    
    % Plot the normalized error (error per one vector)
    fprintf('KSVD():: Iteration %i/%i, |.|_F norm %0.2f\n', i, iters, norm(X - D*A, 'fro'));
    fprintf('Estimated running time left %0.2f minutes \n\n', (iters - i)*time/60);
end

end


% ------------------------------------------------------------------------
% Find a new better dictionary element, i.e. atom vector d
% ------------------------------------------------------------------------
function [d, A] = getd(X,D,j,A)

% Data indices which use the j'th dictionary D atom d_j
indices = find(A(j,:));

% No one uses this atom d_j
if (length(indices) < 1)
    
    % Reconstruction error
    E = X - D*A;
    [~, ind] = max( sum(E.^2) );
    
    % Select a better dictionary atom, normalize, multiply with the sign
    d = X(:,ind);
    d = d / sqrt(d'*d);
    d = d * sign(d(1));
    
    % Update the coefficient matrix
    A(j,:) = 0;
else
    
    % Restrict E_k by choosing only the columns which used d_j
    tmpA = A(:, indices);
    tmpA(j,:) = 0;
    
    % Reconstruction error
    E = X(:,indices) - D*tmpA;
    
    % New dictionary atom is found by using SVD (singular value decomp.)
    [U,D,V] = svd(E, 'econ');
    d = U(:,1);
    
    % We are minimizing: || X - DA ||_F ^ 2. A rank one
    % approximation is done with using the largest singular value.
    
    % Update the coefficient matrix
    A(j,indices) = D(1,1) * V(:,1)';
end

end
% Orthogonal Matching Pursuit (OMP) (greedy ell_0-optimization)
% for sparse signal representations
% ------------------------------------------------------------------------
% 
% Input:   D  =  Dictionary (d x K), with columns normalized
%          X  =  Signal matrix (d x N) to be represented with model X = DA
%          L  =  Maximum number of non-zero coefficients for each signal
%
% Output:  A  =  Sparse gain coefficient matrix (K x N)
% 
% 
% Pati Y, Rezaiifar R, Krishnaprasad P.
% "Orthogonal matching pursuit: recursive function
% approximation with applications to wavelet decomposition",
% Signals, Systems and Computers, 1993.
%
% http://www.isr.umd.edu/~krishna/images/pati_reza_psk.pdf
%
% Mikael Mieskolainen, mikael.mieskolainen@tut.fi, 2011
% ------------------------------------------------------------------------

function A = OMP(D, X, L)

[~, K] = size(D);
[~, N] = size(X);

A = sparse(K, N);

% Go through all N measurement vectors 
for k = 1:N
    
    x = X(:,k);                    % The measurement vector
    
    a = [];
    residual = x;
    ind = zeros(L,1);
    
    for j = 1:L
        
        proj = D'*residual;        % Inner products: <D_column, residual>
        [~, pos] = max(abs(proj)); % Find the index with most correlation
        ind(j) = pos(1);           % this corresponds to column in D
        
        % Solve a least-squares problem to obtain a new estimate
        %a = pinv(D(:,ind(1:j))) * x;
        a = D(:,ind(1:j)) \ x;     % if '\' fails, then use pinv()
        
        % Calculate the new residual
        residual = x - D(:,ind(1:j))*a;
    end
    
    temp = zeros(K,1);             % Create a vector of length K
    temp(ind(1:j)) = a;            % Fill the non-zero coefficients with a
    A(:,k) = sparse(temp);         % Put the estimate to coefficient matrix
end

end
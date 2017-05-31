% CELP Codebook search for optimum vector v and corresponding gain u
% ------------------------------------------------------------------------
%
% Input:    C   =   Codebook matrix
%           H   =   Convolution matrix
%           d   =   Target vector
%
% Output:   v   =   Optimum excitation vector in the codebook C
%          mu   =   Corresponding gain
%         ind   =   Index in the codebook of the best vector
% ------------------------------------------------------------------------
% Target: min ||muHv - d||^2
%
% For adaptive codebook: S = W*H, mu = alpha, v = a
% For fixed codebook:    S = W*H, mu = phi,   v = f
%
% Mikael Mieskolainen, mikael.mieskolainen@tut.fi, 2011
% ------------------------------------------------------------------------

function [v,mu,ind] = codebook_search(C,H,d)

d = d(:);                    % Make sure it is column vector
[N,M] = size(C);             % Codebook dimensions

if all(C(:) == 0)
    v = zeros( 1, M);
    mu = 0;
    ind = 1;
    return;
end

HC = H*C.';                  % Calculate all synthesis vectors
alphas = d.'*HC./sum(HC.^2); % Optimal alphas

% Squared error, end the corresponding best (minimum error) index
err = sum( (d*ones(1,N) - ones(M,1)*alphas.*HC).^2 );
[~, ind] = min(err);

% The best codebook vector and corresponding gain
v = C(ind,:);
mu = alphas(ind);

end
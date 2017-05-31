% Vector quantization based on a codebook, using L2-distance
% ------------------------------------------------------------------------
%
% Input:       cent  :  Codebook (dxK), d is the dimension
%                 v  :  Input vector to be quantized
%
% Output:        vq  :  Closest vector to v in codebook
%               inf  :  The codebook index
%
% Mikael Mieskolainen, mikael.mieskolainen@tut.fi, 2011
% ------------------------------------------------------------------------

function [vq, ind] = VQ_quant(cent, v)

v = v(:);

% Calculate the distances
dist = sum((cent - v*ones( 1, size( cent,2))).^2);

% Find the minimum vector
[~, ind] = min(dist);
vq = cent(:, ind);

end
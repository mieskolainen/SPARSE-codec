% Scalar Quantization
% ------------------------------------------------------------------------
% 
% Input:     CB  =  Codebook (1xN)
%             v  =  Value to be quantized
% Output:     q  =  Quantized value
%           ind  =  Codebook index
% 
% Mikael Mieskolainen, mikael.mieskolainen@tut.fi, 2011
% ------------------------------------------------------------------------

function [q, ind] = VQ1(CB, v)

d = abs(CB - repmat(v, 1, length(CB)));
[~, ind] = min(d);
q = CB(ind);

end
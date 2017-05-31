% Calculate SNR of y and quantized yq
%
%      input:   y = Original signal
%              yq = Quantized signal
%
% Mikael Mieskolainen, mikael.mieskolainen@tut.fi, 2011

function s = SNR(y, yq)

e = y(:)-yq(:);
s = 10*log10( sum(y.^2) / sum(e.^2) );

end
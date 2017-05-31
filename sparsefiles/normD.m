% Normalize dictionary D, column by column with L2-norm to unit length
% ------------------------------------------------------------------------
% 
% Input:       D  =  Dictionary (d x K), where d is dimension, 
%                    K is number of atoms
% Output:      D  =  Normalized dictionary
%
% Mikael Mieskolainen, mikael.mieskolainen@tut.fi, 2011
% ------------------------------------------------------------------------

function D = normD(D)

D = D./repmat(sqrt(sum(D.^2,1)),size(D,1),1);

end
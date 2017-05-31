% Main example for using the 8 kHz sample rate CELP codec
% 
% See:
% M. R. Schroeder and B. S. Atal, "Code-excited linear prediction (CELP), 1985
%
%
% Mikael Mieskolainen, mikael.mieskolainen@tut.fi, 2011

clear;

addpath commonfiles;
addpath celpfiles;

train = false;
%train = true;


%% Train the CELP codec and setup global parameters
% CODEC is hard coded for fs = 8 kHz, but that can be easily changed

if (train)
    param.LSF_bits = 9;   % LSF codebook size
    param.AC_bits = 7;    % Adaptive codebook delays size, max 7
    param.FC_bits = 10;   % Fixed codebook dictionary size
    param.GA_bits = 5;    % Adaptive codebook gains
    param.GF_bits = 5;    % Fixed codebook gains

    % Train the codebooks
    traincelp(param, 'celp_param.mat', './trainwav/');
end

%% Do the encoding and decoding

y = encoder('./testwav/test.wav', 'celptest.bin', 'celp_param.mat');
yq = decoder('celptest.bin', 'celp_param.mat');


%% Here one can truly verify results

fprintf('SNR = %0.3f \n', SNR(y, yq));
figure;
plot(y); 
xlabel('sample', 'interpreter','latex'); 
ylabel('amplitude', 'interpreter','latex');
hold on;
plot(yq, 'r');
legend('Original', 'Decoded');
title('CELP Codec','interpreter','latex');

%soundsc(yq, 8000);


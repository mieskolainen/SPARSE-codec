% Main example for using the SPARSE codec
%
%
% Mikael Mieskolainen, mikael.mieskolainen@tut.fi, 2011
% ------------------------------------------------------------------------
% 
% SPARSE codec is based on
%
% K-SVD to train an overcomplete dictionary based on MDCT
% (Modified Discrete Cosine Transform) vectors
% + OMP (Orthogonal Matching Pursuit) \ell_0 greedy sparse coding
%
% This approach is new and perhaps not used before in time
% domain signal coding, such as speech coding.
%
% NOTE! This algorithmic combination does not assume any properties of
% the signal (such as human speech), so may work with any kind of signals.
% It is just matter of training dictionaries with K-SVD.
% Extensions could include multiple dictionaries and
% multidimensional signals (2-D, N-D, signals on spheres etc.)
% with different basis transforms than MDCT. It is possible also to use
% just raw signals with K-SVD, without any pre-transforms.
%
% ------------------------------------------------------------------------
%
% Every MDCT vector is represented as a linear combination of number of L
% trained atom vectors from the dictionary D. This can be seen as a signal
% adaptive and overcomplete (thus, not unique) linear transform.
% These L vectors and their gains are selected from dictionary D with the
% greedy OMP ell_0 sparse algorithm.
% 
% - Model is:   X  =  DA, where A is the sparse gain matrix, where every
%                         column has L non-zero elements. X is our
%                         measurements vector matrix and D is
%                         the dictionary with atoms as columns.
%
% - MDCT is based on 50 % overlap with sine window, which obeys the TDAC
%  (time domain antialias cancellation) rules.
%
% ------------------------------------------------------------------------
%
% NOTE! MDCT transform (files in ./mdctlib folder) are libraries from:
%                                  http://www.ee.columbia.edu/~marios/
%
% Unzip the mdctlib.zip to ./mdctlib folder.
% Error with 'linframe' etc. will occur without this.
% ------------------------------------------------------------------------

clear;

addpath commonfiles;
addpath sparsefiles;
addpath mdctlib;

train = false;
%train = true;


%% Train and setup global parameters, Note! This is slow.
% Took about few hours on 4-core 2.4 GHz Q6600 with every
% core in use. Needs to be done only once, though.

if (train)
    param.L = 8;          % Sparseness \ell_0, affects quality vs. compress ratio
    param.CB_bits =  10;  % Dictionary bits, dictionary size K is 2^CB_bits
    param.G_bits = 6;     % Bits per gain coefficient
    param.WIN_SIZE = 256; % Window length
    param.Fs = 8000;      % Sampling frequency

    iterations = 25;      % Iterations in K-SVD, should be at least > 25
    samples = 20000;      % Number of training vectors in K-SVD, should be at least
                          % thousands or tens of thousands

    trainsparse(param, './trainwav/', 'sparse_param.mat', iterations, samples);
end


%% Encode and decode, this is fast.

y = encodesparse('./testwav/test.wav', 'sparsetest.bin', 'sparse_param.mat');
yq = decodesparse('sparsetest.bin', 'sparse_param.mat');


%% Here one can once more verify the results

fprintf('SNR = %0.3f \n', SNR(y, yq));

figure;
plot(y); hold on;
plot(yq, 'r');
xlabel('sample','interpreter','latex');
ylabel('amplitude','interpreter','latex');
legend('Original', 'Decoded');
title('SPARSE codec', 'interpreter', 'latex');

%soundsc(yq, 8000);


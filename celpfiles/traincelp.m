% Train the CELP codec (fixed with Fs = 8000 Hz)
% ------------------------------------------------------------------------
% 
% Input: param   =   Parameter struct, for example:
%                    param.LSF_bits = 8;, LSF dictionary size
%                    param.AC_bits = 7; , Adaptive codebook delays size, max 7
%                    param.FC_bits = 10;, Fixed codebook dictionary size
%                    param.GA_bits = 4; , Adaptive codebook gains
%                    param.GF_bits = 4; , Fixed codebook gains
%
%        output_file = Filename to save the parameters
%        trainpath   = Path to .wav trainfiles, e.g. 'C:/Trainfiles/' 
%
% Mikael Mieskolainen, mikael.mieskolainen@tut.fi, 2011
% ------------------------------------------------------------------------

function traincelp(param, output_file, trainpath)

LSF_bits = param.LSF_bits;
AC_bits = param.AC_bits;
FC_bits = param.FC_bits;
GA_bits = param.GA_bits;
GF_bits = param.GF_bits;

%------------------------------------------------------------------------
% LPC order fixed
P = 10;

% ------------------------------------------------------------------------
% Get all the train files to struct array
trainfiles = dir([trainpath '/*.wav']);
if (isempty(trainfiles))
    error('Trainpath does not include any wav file!'); 
end

fprintf('Training CELP codec with %d .wav files from %s\n', ...
        length(trainfiles), trainpath);

% ------------------------------------------------------------------------

% Construct the fixed codebook (stochastic codebook)
FC = randn(2^FC_bits, 40);

% Normalize to unit energy
for i = 1:size( FC,1)
    FC( i,:) = FC(i,:)/sum( FC(i,:).^2);
end

% ------------------------------------------------------------------------

% Construct the indices for the adaptive codebook
% Delays in taps [16 ... 144] ~ [55.5 Hz ... 500 Hz] in frequencies
ACind = zeros(2^AC_bits, 40);
for i = 1:length(ACind)
    delay = 15 + i;
    delayvec = -delay + mod(0:39, delay); % Delays relative to this sample  
    ACind(i,:) = delayvec;
end

% ------------------------------------------------------------------------
% Codebook of LSF vectors, times 5 because of split vector quantization SVQ
LSF_cents = train_LSF(trainpath, trainfiles, P, LSF_bits);

% AT this point we will use these temporary codebooks, these will be replaced
% with K-means trained at lines 85 ->
AG = linspace(-2, 2, 1024);
FG = linspace(-2, 2, 1024);

save(output_file,'FC_bits','AC_bits','LSF_bits','GA_bits','GF_bits','P','FC',...
                 'ACind','LSF_cents','AG','FG')

% ------------------------------------------------------------------------
% Finally train the gain codebooks with K-means 
all_alfas = []; all_phis = [];

% Encode the wave files
for i = 1:length(trainfiles)
    [~, alfas, phis] = encoder([trainpath '/' trainfiles(i).name],...
                                       'tmp.bin', output_file);
    all_alfas = [all_alfas alfas];
    all_phis = [all_phis phis];
end

% Train gain codebooks now with K-means
[AG, ~, ~] = K_means(all_alfas, 2^GA_bits);
[FG, ~, ~] = K_means(all_phis, 2^GF_bits);

% Save again
save(output_file,'FC_bits','AC_bits','LSF_bits','GA_bits','GF_bits','P','FC',...
                 'ACind','LSF_cents','AG','FG')

end

% ------------------------------------------------------------------------
% Train LSF codebook (Fixed with Fs = 8000 Hz)
% ------------------------------------------------------------------------
%
% Input:    trainpath   =  Path to trainfiles
%           trainfiles  =  Struct array of trainfilenames
%                    p  =  Order of LPC analysis
%                 bits  =  Number of bits, size of codebook will be 2^bits
%
% Output:            C  =  LSF codebook (P x 2^bits)
%
% Mikael Mieskolainen, mikael.mieskolainen@tut.fi, 2011
% ------------------------------------------------------------------------

function cents = train_LSF(trainpath, trainfiles, p, bits)

% Get all the train files
signal = [];
for i = 1:length(trainfiles)
    [this_signal, Fs] = audioread([trainpath '/' trainfiles(i).name]);
    if (Fs ~= 8000)
        error('Wave file should be with Fs = 8000 Hz');
    end
    signal = [signal; this_signal];
end

win_length = 20e-3 / (1/8000);     % 20 ms window
win_jump = 2.5e-3 / (1/8000);      % 2.5 ms jump between windows
window = hann(win_length);         % Hanning window

% Indexes of the signal
indexes = 1:win_jump:length(signal) - win_length;

% Reflection coefficients matrix
reflmat = zeros(p, length(indexes));

% Go through the signal in frames of 20 ms and window with a Hanning
i = 1;
for n = indexes
    
    % Windowed frame
    frame = signal(n:n+win_length-1) .* window;
    
    % Calculate LPC's of order p
    a = lpc(frame, p);
    
    % Convert to LSF coefficients and save to matrix
    reflmat(:,i) = poly2lsf(a);
    i = i + 1;
end

% Now we have the trainingvectors (reflmat)
% Train 5 codebooks of dimension 2, each with K vectors
% If we use e.q bits = 7, number of centroids K = 128
K = 2^bits;
cents = cell(5,1);

% Run the K-means for the first 2 elements for the matrix generated above
for i = 1:5
    [cents{i}, ~, ~] = K_means(reflmat(2*i-1:2*i,:), K);
end

end
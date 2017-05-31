% Train the SPARSE codec
% ------------------------------------------------------------------------
%
% Input:        param   =  Parameters struct, for example:
%
%  param.L = 8;          , Sparseness \ell_0, affects quality vs. compress ratio
%  param.CB_bits = 10;   , Dictionary bits, dictionary size K is 2^CB_bits
%  param.G_bits = 6;     , Gain bits
%  param.WIN_SIZE = 256; , Window length
%  param.Fs = 8000;      , Sampling frequency
%
%            trainpath  =  Path to train wave files
%           iterations  =  Number of iterations in K-SVD
%             samples   =  Number of sample vectors in K-SVD training
%
% Mikael Mieskolainen, mikael.mieskolainen@tut.fi, 2011
% ------------------------------------------------------------------------

function trainsparse(param, trainpath, param_file, iterations, samples)

% Window size and sampling frequency
WIN_SIZE = param.WIN_SIZE;
Fs = param.Fs;

% Codec parameters
CB_bits = param.CB_bits;    % Codebook bits
G_bits = param.G_bits;      % Gain bits
K = 2^CB_bits;              % -> Codebook size
L = param.L;                % Sparseness

% Get train vectors
files = dir([trainpath '/*.wav']);
filelist = cell(length(files), 1);
for i = 1:length(files)
    filelist{i} = [trainpath '/' files(i).name];
end

% Get speech sample vectors
X = getMDCT(filelist, WIN_SIZE, Fs, samples);

% Train dictionary D with K-SVD
[D, A] = ksvd(X, K, iterations, 'Gaussian', L);

% GCB (Gain CodeBook), uniform quantization on interval [-x ... x].
% Uniform quantization is done faster way in a real system, but here we
% just use the 1D-codebook approach.
%GCB = linspace(-0.7, 0.7, 2^G_bits);

% Train gain codebook now with K-means
fprintf('Training Gain Codebook with K-means \n');
[GCB, ~, ~] = K_means(full(A(A~=0))', 2^G_bits);

% ------------------------------------------------------------------------
% These are frequencies/probabilities for arithmetic entropy coding

% Train the index frequency (probability) table
fprintf('Training dictionary index probabilities \n');
freq_i = zeros(size(D,2), 1);
for i = 1:size(A,1)
    freq_i(i) = length(find(A(i,:))); % Find the number of non-zero weights
end

% Train the weight frequency (probability) table
fprintf('Training dictionary atom weight probabilities \n');
freq_w = zeros(length(GCB), 1);
for i = 1:size(A,1)
    for j = 1:size(A,2)
        if (A(i,j) ~= 0)
            [~,ind] = VQ1(GCB, A(i,j) );
            freq_w(ind) = freq_w(ind) + 1;
        end
    end
end

% ------------------------------------------------------------------------

fprintf('Training done! \n');

% Save all necessary
save(param_file, 'D', 'GCB', 'CB_bits', 'G_bits', 'L', ...
                 'WIN_SIZE', 'Fs', 'freq_w', 'freq_i');
end


% ------------------------------------------------------------------------
% Collect Modified Discrete Cosine Transform (MDCT) vectors
% ------------------------------------------------------------------------
% 
% Input:   filelist  =  List of files to be used in training (cell)
%        win_size_s  =  Window length in samples
%                fs  =  Sampling frequency
%                 M  =  Number of sample vectors
%                 
% Output:         X  =  Measurements matrix (d x M)
%
% Mikael Mieskolainen, matti.m.mieskolainen@tut.fi
% ------------------------------------------------------------------------

function allX = getMDCT(filelist, win_size_s, fs, M)

% Here we put all the training vectors
allX = zeros(win_size_s/2, M);

% Collect all the sample vectors
sample_id = 1;
sample_vectors_per_file = floor(M / length(filelist)) + 1;

for file_id = 1:length(filelist)
    
    % Read files randomly again and again, slow but is enough here
    [y, fs_this] = audioread( filelist{file_id} );
    if (fs_this ~= fs)
        erros(['Wrong Fs with file' filelist{file_id}]);
    end
    
    % Select random locations
    for i = 1:sample_vectors_per_file
        r = floor(rand(1)*(length(y) - win_size_s + 1)) + 1;
        ynow = y(r:r + win_size_s - 1);

        fx = winit(ynow,'sinewin');            % Window the signal
        if (sample_id <= M)
            allX(:, sample_id) = mdct4(fx);    % Put to the matrix
            sample_id = sample_id + 1;
        end
    end
end

% Normalize all the measurement vectors to have L2-unit norm
% allX = normD(allX);
% NOTE!, THIS NORMALIZATION CAUSES NOISE TO SILENT PARTS (SO DO NOT USE!)

fprintf('Reading of %i MDCT training vectors done \n', M);

end
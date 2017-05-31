% SPARSE codec decoder
% ------------------------------------------------------------------------
%
% Input:     input_file  =  Input binary file
%            param_file  =  Parameter .mat file
%            verbose     =  Print output (true), no output (false)
% 
% Output:            yq  =  Decoded signal vector
% 
% Mikael Mieskolainen, mikael.mieskolainen@tut.fi, 2011
% ------------------------------------------------------------------------

function yq = decodesparse(input_file, param_file, verbose)

if (nargin < 3)
   verbose = true; 
end

% Load the coded parameters
load(param_file);

% Load the binary file
fid = fopen(input_file, 'r');

A_size = fread(fid, [2 1], 'uint32');                                % Signal size (header)
pad_length = fread(fid, 1, 'uint32');                                % Pad length (header)

%positions = fread(fid, [A_size(2)*L 1], ['ubit' int2str(CB_bits)]); % Dictionary positions
%gains = fread(fid, [A_size(2)*L 1], ['ubit' int2str(G_bits)]);      % Gains

index_code_length = fread(fid, 1, 'uint32');                         % Length of index bit vector
weight_code_length = fread(fid, 1, 'uint32');                        % Length of weight bit vector

index_L = fread(fid, 1, 'uint32');                                   % Number of positions
index_code = fread(fid, [index_code_length 1], 'ubit1');             % Dictionary positions
weight_code = fread(fid, [weight_code_length 1], 'ubit1');           % Gains

fpad = fread(fid, [pad_length 1], 'float');                          % Remaining file part

% Fix MATLAB indexing
%positions = positions + 1;
%gains = gains + 1;

% Close the file
fclose(fid);

% ------------------------------------------------------------------------
% Arithmetic Decoding

positions = arithdeco(index_code, freq_i, index_L);
gains = arithdeco(weight_code, freq_w, index_L);


% ------------------------------------------------------------------------
tic;

% Reconstruct matrix A
A = reconstruct_A(positions, gains, L, GCB, A_size);

% Inverse MDCT
fy = imdct4(D*A);
fy = winit(fy,'sinewin');               % Rewindow
yq  = linunframe(fy, WIN_SIZE/2, fpad); % Overlap add (OLA) to get the decoded signal

decode_time = toc;

if (verbose)
    fprintf('\n----------------------------------------------------------------\n');
    fprintf('SPARSE CODEC - Decoding results (%s -> yq) \n', input_file);
    fprintf('----------------------------------------------------------------\n');
    fprintf('File length: %0.1f s / %d samples \n', length(yq)/Fs, length(yq));
    fprintf('Decoding time: %0.2f s, %0.2f x realtime \n', decode_time, (length(yq)/8000) / decode_time);
    fprintf('----------------------------------------------------------------\n\n');
end

end


% ------------------------------------------------------------------------
% Reconstruct gain matrix A
% -----------------------------------------------------------------------
function A = reconstruct_A(positions, gains, L, GCB, A_size)

A = zeros(A_size');
col = 1;

for i = 1:L:length(positions)
    
    a = positions(i:L+i-1); % Sparse column of original A
    g = gains(i:L+i-1);     % Gain codebook indices
    
    for j = 1:length(a)
        A(a(j), col) = GCB(g(j)); % Select the gain from codebook 
    end
    col = col + 1;
end

end
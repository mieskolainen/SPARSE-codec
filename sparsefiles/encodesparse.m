% SPARSE codec encoder
% ------------------------------------------------------------------------
% 
% Input:     input_file   =   Input wave filename
%           output_file   =   Output binary filename
%            param_file   =   Parameters .mat file
%
% Output:             y   =   Original wave signal as vector
%
% Mikael Mieskolainen, mikael.mieskolainen@tut.fi, 2011
% ------------------------------------------------------------------------

function y = encodesparse(input_file, output_file, param_file)

% Load the coded parameters
load(param_file);

% Read a speech file
[y, Fs_this] = audioread(input_file);
if (Fs_this ~= Fs)
   error(['Wave file is not with Fs = ' int2str(Fs)]); 
end

tic;

% Window the signal, 50 % overlap
[fx, fpad] = linframe(y, WIN_SIZE/2, WIN_SIZE, 'sym');  
fx = winit(fx,'sinewin'); % Sinewindow obeys the TDAC 
                          % (time domain antialias cancellation) conditions

% Do the Modified Discrete Cosine Transform (MDCT)
FX = mdct4(fx);

% Encode using the greedy L_0 Orthogonal Matching Pursuit
A = OMP(D, FX, L);

% ------------------------------------------------------------------------

% Now do the bit-level coding
positions = zeros(L*size(A,2), 1);
gains = zeros(L*size(A,2), 1);

i = 1;
for n = 1:size(A,2)
    for m = 1:size(A,1)
        if (A(m,n) ~= 0)
            positions(i) = m;                 % The column in dictionary D
            [~, gains(i)] = VQ1(GCB, A(m,n)); % 1D-vector quantization
            i = i + 1;
        end
    end
end

% Fix the MATLAB problematic indexing
%positions = positions - 1;
%gains = gains - 1;

% ------------------------------------------------------------------------
% Arithmetic Encoding

index_code = arithenco(positions, freq_i);
weight_code = arithenco(gains, freq_w); 

fprintf('Index: Length with Rissanen %d, without %d \n', length(index_code), CB_bits*length(positions));
fprintf('Weights: Length with Rissanen %d, without %d \n', length(weight_code), G_bits*length(gains));

% Huffman Entropy Encoding
%{
set(0,'RecursionLimit', 2000); % Set recursion limit larger
[pos_dict,avglen] = huffmandict(1:2^CB_bits, freq_i / sum(freq_i) );
[weight_dict,avglen] = huffmandict(1:2^G_bits, freq_w / sum(freq_w) );

index_code_huff  = huffmanenco(positions, pos_dict);
weight_code_huff = huffmanenco(gains, weight_dict);

fprintf('Index: Length with Huffman %d, without %d \n', length(index_code_huff), CB_bits*length(positions));
fprintf('Weights: Length with Huffman %d, without %d \n', length(weight_code_huff), G_bits*length(gains));
%}

% ------------------------------------------------------------------------

% Save the binary
fid = fopen(output_file, 'w');
fwrite(fid, size(A), 'uint32');                     % Signal size (header)
fwrite(fid, length(fpad), 'uint32');                % Pad length (header)

%fwrite(fid, positions, ['ubit' int2str(CB_bits)]); % Dictionary positions
%fwrite(fid, gains, ['ubit' int2str(G_bits)]);      % Gains

fwrite(fid, length(index_code), 'uint32');          % Dictionary positions
fwrite(fid, length(weight_code), 'uint32');         % Gains

fwrite(fid, length(positions), 'uint32');           % Number of positions
fwrite(fid, index_code, 'ubit1');                   % Dictionary positions bits
fwrite(fid, weight_code, 'ubit1');                  % Gains bits

fwrite(fid, fpad, 'float');                         % Remaining file part

fclose(fid);
% ------------------------------------------------------------------------

% Get the file properties (bytes size etc.)
comp = dir(output_file);
orig = dir(input_file);
encode_time = toc;

% Do the decoding to get SNR
yq = decodesparse(output_file, param_file, false);

fprintf('\n----------------------------------------------------------------\n');
fprintf('SPARSE CODEC - Encoding results (%s -> %s) \n', input_file, output_file);
fprintf('----------------------------------------------------------------\n');
fprintf('Compression ratio: %0.3f \n', comp.bytes / orig.bytes);
fprintf('Bit rate: %0.3f kbit/s \n', comp.bytes*8/(length(y)/Fs)/1024);
fprintf('SNR: %0.1f dB \n', SNR(y, yq));
fprintf('File length: %0.1f s / %d samples \n', length(y)/Fs, length(y));
fprintf('Encoding time: %0.2f s, %0.2f x realtime \n', encode_time, (length(y)/Fs) / encode_time);
fprintf('----------------------------------------------------------------\n\n');

end
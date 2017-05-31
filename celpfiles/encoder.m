% CELP encoder (Fixed for Fs = 8000 Hz)
% ------------------------------------------------------------------------
% 
% Input:   input_file  =  Wav file full path, e.g. 'test.wav'
%         output_file  =  Output binary file full path, e.g 'output.bin'
%     parameters_file  =  Parameters .mat file path 
% 
% Output:           y  =  Original wave file vector
%               alfas  =  Adaptive Codebook gains (for training purposes)
%                phis  =  Fixed Codebook gains (for training purposes)
%
% Mikael Mieskolainen, mikael.mieskolainen@tut.fi, 2011
% ------------------------------------------------------------------------

function [y, alfas, phis] = encoder(input_file, output_file, parameters_file)

% Load parameters
load(parameters_file);

% Load the wave file
[y, Fs] = audioread(input_file);
if (Fs ~= 8000)
   error('Wav file should be with Fs = 8000 Hz'); 
end

y = y(:);
y_length = length(y);

% ------------------------------------------------------------------------
% Create analysis window
awin = hanning(400); awin = awin(1:200);      % Long half
tmpwin = hanning(80);                         % Temp window
awin(201:240) = tmpwin(41:80);                % Final analysis window

% ------------------------------------------------------------------------

subf0 = zeros(1,40);               % One subframe of zeros
wprev = poly2lsf([1 zeros(1,P)]);  % initial LFS's from initial LPC
Zi = zeros(P, 1);                  % Initial filter memory

n = 1;                             % Running index for main sample index
sfind = 1;                         % Running index for each subframe
lindex = 1;                        % Running index for each split VQ

% ------------------------------------------------------------------------

yq = zeros(1, y_length);           % quantized signal goes here
yr = zeros(1, y_length);           % we construct the residual here

% ------------------------------------------------------------------------
% These we will ENCODE

LSFinds = ones(ceil(length(y)/160) * 5, 1);    % LSF indices
delays = ones(ceil(length(y)/40), 1);          % AC delays
fcinds = ones(ceil(length(y)/40), 1);          % FC indices
alfainds = ones(ceil(length(y)/40), 1);        % AC Gains
phinds   = ones(ceil(length(y)/40), 1);        % FC Gains

% ------------------------------------------------------------------------

tic;

while ( n + 240 <= y_length)              % loop through signal
    
    af = y( n+(0:length(awin)-1) ).*awin; % analysis window
    
    a = lpc(af, P);                       % LPC analysis
    w = poly2lsf( a);                     % ...to lsf's
    
    % Split Vector Quantization of LSF coefficients
    % Quantize sub vectors: {1,2},{3,4},{5,6},{7,8},{9,10}
    for i = 1:5
        [w(i*2-1:i*2), LSFinds(lindex)] = VQ_quant(LSF_cents{i}, w(i*2-1:i*2));
        lindex = lindex + 1;
    end
    
    % Go through all the subframes
    for subind = 1:4
        
        substart = n + subind*40;        % First sample of the subframe
        
        % Get LPC polynomial for subframe using linear interpolation
        lsfi = (1-subind/4)*wprev+(subind/4)*w; % LSF weighting for subframe
        ai = lsf2poly(lsfi); % this is the LPC polynomial for this subframe
        
        % First get the zero input response (ZIR)
        zir = filter(1, ai, subf0, Zi);
        
        % Now get target vector for this subframe
        d = y(substart + (0:39)).';
        
        % Remove the zero input response (ZIR)
        d = d - zir; % this is the ideal output of the synthesis filter
                     % without memory
        
        % Calculate the convolution matrix
        h = impz(1, ai, 40);  % Impulse response of the synthesis filter
        H = convmtx(h, 40);
        H = H(1:40,1:40);     % Chop the band matrix to approximate the IIR
        
        % Construct the adaptive codebook from the residual
        if (n > 144)
            AC = yr(substart + ACind); % Index using matrices
        else
            AC = zeros(2^AC_bits, 40); % For the first subframe, use zeros
        end
        
        % Adaptive codebook search and gain quantization and sign bit
        [av, alfa, delays(sfind)] = codebook_search(AC, H, d);
        alfas(sfind) = alfa; % (Save for training function before quantization)
        [alfa, alfainds(sfind)] = VQ1(AG, alfa);
        
        % Update the target vector
        d = d - (alfa*H*av(:)).';
        
        % Fixed codebook search and gain quantization and sign bit
        [fv, phi, fcinds(sfind)] = codebook_search(FC, H, d);
        phis(sfind)  = phi; % (Save for training function before quantization)
        [phi, phinds(sfind)] = VQ1(FG, phi);
        
        % ----------------------------------------------------------------
        
        % The residual signal for this subframe is alfa*av + phi*fv
        res = alfa*av + phi*fv; 
        yr(substart + (0:39)) = res; % update the residual
        
        % The output signal is now res filtered through the synthesis 
        % filter. Update the filter memory at the same time.
        [subout, Zi] = filter(1, ai, res, Zi);
        yq(substart + (0:39)) = subout;
        
        sfind = sfind + 1;
    end
    
    wprev = w; % Update previous lsf's
    n = n + 160;
end

% Before saving, fix the MATLAB problematic indexing to avoid overflow
LSFinds = LSFinds - 1;
delays = delays - 1;
fcinds = fcinds - 1;
alfainds = alfainds - 1;
phinds = phinds - 1;

% Save the binary
fid = fopen(output_file, 'w');

fwrite(fid, y_length, 'uint64');                  % Signal length
fwrite(fid, LSFinds, ['ubit' int2str(LSF_bits)]); % LSF VQ indices
fwrite(fid, delays, ['ubit' int2str(AC_bits)]);   % Adaptive codebook (AC) delays
fwrite(fid, fcinds, ['ubit' int2str(FC_bits)]);   % Fixed codebook (FC) indices
fwrite(fid, alfainds, ['ubit' int2str(GA_bits)]); % Fixed codebook (FC) gain indices
fwrite(fid, phinds, ['ubit' int2str(GF_bits)]);   % Adaptive codebook (AC) gain indices

fclose(fid);

% Get the file properties (bytes size etc.)
comp = dir(output_file);
orig = dir(input_file);
encode_time = toc;

fprintf('\n----------------------------------------------------------------\n');
fprintf('CELP CODEC - Encoding results (%s -> %s) \n', input_file, output_file);
fprintf('----------------------------------------------------------------\n');
fprintf('Compression ratio: %0.3f \n', comp.bytes / orig.bytes);
fprintf('Bit rate: %0.3f kbit/s \n', comp.bytes*8/(y_length/8000)/1024);
fprintf('SNR: %0.1f dB \n', SNR(y, yq));
fprintf('File length: %0.1f s / %d samples \n', y_length/8000, y_length);
fprintf('Encoding time: %0.1f s, %0.2f x realtime \n', encode_time, (y_length/8000) / encode_time);
fprintf('----------------------------------------------------------------\n\n');

% %%
% figure
% hist( alfas, 50);
% title('histrogram of adaptive codebook gains');
% 
% figure
% hist( phis, 50);
% title('histrogram of fixed codebook gains');

end
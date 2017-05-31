% CELP decoder (Fixed for Fs = 8000 Hz)
% ------------------------------------------------------------------------
% 
% Input:       input_file  =  Binary file full path, e.g 'test.bin'
%         parameters_file  =  Parameters .mat file path 
% 
% Output:              yq  =  Decoded speech vector
%
% Mikael Mieskolainen, mikael.mieskolainen@tut.fi, 2011
% ------------------------------------------------------------------------

function yq = decoder(input_file, parameters_file)

% Load the parameters
load(parameters_file);

% Load the binary file
fid = fopen(input_file, 'r');

y_length = fread(fid, 1, 'uint64');                                        % Signal length
LSFinds = fread(fid, [ceil(y_length/160)*5 1], ['ubit' int2str(LSF_bits)]);% LSF VQ indices
delays = fread(fid, [ceil(y_length/40) 1], ['ubit' int2str(AC_bits)]);     % Adaptive codebook (AC) delays
fcinds = fread(fid, [ceil(y_length/40) 1], ['ubit' int2str(FC_bits)]);     % Fixed codebook (FC) indices
alfainds = fread(fid, [ceil(y_length/40) 1], ['ubit' int2str(GA_bits)]);   % Fixed codebook (FC) gain indices
phinds = fread(fid, [ceil(y_length/40) 1], ['ubit' int2str(GF_bits)]);     % Adaptive codebook (AC) gain indices

fclose(fid);

% Before continuing, fix the MATLAB indexing 0 -> 1
LSFinds = LSFinds + 1;
delays = delays + 1;
fcinds = fcinds + 1;
alfainds = alfainds + 1;
phinds = phinds + 1;

% ------------------------------------------------------------------------

w = zeros(P,1);                    % LSF coefficients
wprev = poly2lsf([1 zeros(1,P)]);  % initial LFS's from initial LPC
Zi = zeros(P, 1);                  % Initial filter memory

n = 1;                             % Running index for main sample index
sfind = 1;                         % Running index for each subframe
lindex = 1;                        % Running index for each split VQ

% ------------------------------------------------------------------------

yq = zeros(1, y_length);           % quantized signal goes here
yr = zeros(1, y_length);           % we construct the residual here

% ------------------------------------------------------------------------

tic;

while (n + 240 <= y_length)  % loop through signal
    
    % Construct Split Vector Quantization of LSF coefficients
    % Quantize sub vectors: {1,2},{3,4},{5,6},{7,8},{9,10}
    for i = 1:5
        w(i*2-1:i*2) = LSF_cents{i}(:,LSFinds(lindex));
        lindex = lindex + 1;
    end
    
    % Go through all the subframes
    for subind = 1:4
        
        substart = n + subind*40; % first sample of the subframe
        
        % Get LPC polynomial for subframe using linear interpolation
        lsfi = (1-subind/4)*wprev + (subind/4)*w; % lsf weighting for subframe
        ai = lsf2poly(lsfi); % this is the LPC polynomial for this subframe
        
        % Construct the adaptive codebook from the residual
        if (n > 144)
            AC = yr(substart + ACind); % Index nicely in matrix form
        else
            AC = zeros(2^AC_bits, 40); % For the first subframe, use just zeros
        end
        
        % AC vector and gain
        av = AC(delays(sfind),:);
        alfa = AG(alfainds(sfind));
        
        % FC vector and gain
        fv = FC(fcinds(sfind),:);
        phi = FG(phinds(sfind));
        
        % The residual signal for this subframe is alfa*av + phi*fv
        res = alfa*av + phi*fv; 
        yr(substart + (0:39)) = res; % update the residual
        
        % The output signal is residual filtered through the synthesis 
        % filter. Update also the filter memory.
        [subout, Zi] = filter(1, ai, res, Zi);
        yq(substart + (0:39)) = subout;
        
        sfind = sfind + 1;
    end
    
    wprev = w; % Update previous lsf's
    n = n + 160;
end

decode_time = toc;

fprintf('\n----------------------------------------------------------------\n');
fprintf('CELP CODEC - Decoding results (%s -> yq) \n', input_file);
fprintf('----------------------------------------------------------------\n');
fprintf('File length: %0.1f s / %d samples \n', y_length/8000, y_length);
fprintf('Decoding time: %0.1f s, %0.2f x realtime \n', decode_time, (y_length/8000) / decode_time);
fprintf('----------------------------------------------------------------\n\n');

end
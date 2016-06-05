% spectralnoisegating.m - Filters noise in a signal based on correlation
% subtraction as described in the following paper.
%
% Please cite the following paper if you use this code:
% J. M. Inouye, S. S. Blemker, and D. I. Inouye. "Towards Undistorted and 
%    Noise-free Speech in an MRI Scanner: Correlation Subtraction 
%    Followed by Spectral Noise Gating", Journal of the Acoustical Society 
%    of America, 135(3):1019-1022, 2014.
%
% Copyright (C) 2013 Joshua M. Inouye, David I. Inouye
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License version 3 as 
% published by the Free Software Foundation.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
function sngSignal = spectralnoisegating(signal, noiseSample, Fs, noiseReductionCoeff)
%SPECTRALNOISEGATING Cancel noise using spectral noise gating
% Filter a signal using spectral noise gating implemented in SoX (Sound 
% Exchange - <http://sox.sourceforge.net/>).
%
% NOTE: Requires SoX executables.
%
% sngSignal = spectralnoisegating(signal, Fs) returns filtered signal 
% and will use approximately the first second as a noise sample and the
% noise reduction coefficient will be set to 0.5.
%
% sngSignal = spectralnoisegating(signal, Fs, noiseSample) same as above
% except the explicit noiseSample will be used.
%
% sngSignal = spectralnoisegating(signal, Fs, noiseSample, noiseReductionCoeff) 
% same as above except the explicit noiseReductionCoeff will be used.

%% Extract noise sample if none given
if(isempty(noiseSample))
    % Delay noise sample to avoid initialization frequency differences
    delay = floor(Fs/100);
    noiseSample = signal(delay:Fs);
end
if(nargin < 4)
    noiseReductionCoeff = 0.5;
end

%% Write signals to files to be used by SoX
normFactor = max(max(abs(signal)), max(abs(noiseSample)) ) / 0.8;  

noiseProfileFile = [tempdir 'noiseprofile.wav'];
speechSignalFile = [tempdir 'speechsignal.wav'];
wavfilewrite(noiseProfileFile, noiseSample / normFactor, Fs)
wavfilewrite(speechSignalFile, signal / normFactor, Fs)

if(ispc == 1)
    soxExec = 'sox\sox-14.4.1-win32\sox.exe';
elseif(ismac == 1)
    soxExec = 'sox/sox-14.4.1-macosx/sox';
else % Linux or others
    [status, soxExec] = system('which sox');
    soxExec = strtrim(soxExec);
    if(status == 0)
        %fprintf('Using %s as the SoX executable\n', soxExec);
    else
        error(['Cannot find SoX executable.  Have you installed SoX '...
            '<a href="sox.sourceforge.net">(SoX website)</a>?']);
    end
end

% Compute noise profile with SoX
speechProfileFile = [tempdir 'speech.noiseprofile'];
system(sprintf('%s %s -n noiseprof %s',soxExec, noiseProfileFile, speechProfileFile));

% Reduce noise in signal by using noise profile in SoX
cleanSpeechFile = [tempdir 'cleanspeech.wav'];
system(sprintf('%s %s %s noisered %s %f', ...
    soxExec, speechSignalFile, cleanSpeechFile, speechProfileFile, noiseReductionCoeff));

% Read SoX processed wav file and scale to original level
sngSignal = wavfileread(cleanSpeechFile) * normFactor;

% Delete temporary files
delete(noiseProfileFile);
delete(speechSignalFile);
delete(cleanSpeechFile);
delete(speechProfileFile);

%% Helper file write functions for different versions of MATLAB
function wavfilewrite(filename, Y, Fs)
    if(exist('audiowrite')) 
        % Newer versions of MATLAB
        audiowrite(filename, Y, Fs)
    elseif(exist('wavwrite'))
        % Older versions of MATLAB
        wavwrite(Y, Fs, filename);
    else
        error([ 'Cannot write *.wav files because "audiowrite" and "wavwrite" '...
                'do not exist.  Please check your MATLAB installation.'] );
    end
end

function [signal] = wavfileread(filename)
    if(exist('audioread')) 
        % Newer versions of MATLAB
        signal = audioread(filename);
    elseif(exist('wavread'))
        % Older versions of MATLAB
        signal = wavread(filename);
    else
        error([ 'Cannot read *.wav files because "audioread" and "wavread" '...
                'do not exist.  Please check your MATLAB installation.'] );
    end   
end

end
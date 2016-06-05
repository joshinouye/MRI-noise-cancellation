% computenoisereduction.m - Estimates the noise reduction.
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
function [noiseSupp, noiseOnlyNoiseSupp] = computenoisereduction(methodTitle, speechAndNoise, filteredSpeechAndNoise, noiseOnly, filteredNoiseOnly, showNoiseOnly)
%COMPUTENOISEREDUCTION Computes noise reduction and increase in SNR.
% Computes the noise reduction (in dB) and increase in the signal-to-noise
% ratio (SNR) and prints the computations to the terminal.
%
% These calculations are from the following paper:
% E. Bresch, J. Nielsen, K. Nayak, and S. Narayanan, "Synchronized and 
%   noise-robust audio recordings during realtime magnetic resonance 
%   imaging scans," J. Acoust. Soc. Am., vol. 120, no. 4, pp. 1791-1794, 
%   Oct. 2006.
%
% Outputs:
% noiseSupp   - *estimate* of SNR increase from before to after filtering
% noiseOnlyNoiseSupp   - noise suppression (in dB) without voice
%
% Inputs:
% methodTitle            - title of method used (for printout)
% speechAndNoise         - voice and noise
% filteredSpeechAndNoise - filtered signal
% noiseOnly              - (optional) noise only 
% filteredNoiseOnly      - (optional) filtered noise only
%
% [noiseSupp, noiseOnlyNoiseSupp] = computenoisereduction(methodTitle, speechAndNoise, filteredSpeechAndNoise, noiseOnly, filteredNoiseOnly)


%% Calculate noise suppression (noise only)
noiseOnlyNoiseSupp = 10*log10(power(noiseOnly) / power(filteredNoiseOnly));

%% Estimate SNRs to calculate noise suppresion in noisy speech signals
speechPower = power(speechAndNoise) - power(noiseOnly);
snrBeforeFilter = 10*log10( speechPower / power(noiseOnly) );

filteredSpeechPower = power(filteredSpeechAndNoise) - power(filteredNoiseOnly);
snrAfterFilter = 10*log10( filteredSpeechPower / power(filteredNoiseOnly) );

noiseSupp = snrAfterFilter - snrBeforeFilter;
fprintf('%s noise suppression: %.2f dB\n', methodTitle, noiseSupp);

%% Calculate power of signal (subfunction)
function [ power ] = power(signal)
    power = (norm(signal))^2;
end

end
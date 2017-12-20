% correlationsubtraction.m - Filters noise in a signal based on correlation
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
function csSignal = correlationsubtraction(signal, noise)
%CORRELATIONSUBTRACTION Cancel noise using correlation subtraction in time domain
% Cancels background noise by aligning the signal with a similar noise 
% signal using cross correlation and then then subtracting the two waveform
%
% Outputs:
% csSignal   - Filtered signal
%
% Inputs:
% signal          - main signal
% noise           - a noise signal the same length as the main signal
%
% csSignal = computenoisereduction(signal, noise)

%% Subtract noise from speech signal (correlation-subtraction, CS)
% Cross-correlate noise with speech signal
[C,lags]=xcorr(noise, signal);

% Find index of max of cross-correlation sequence
[~,I]=max(C);

% Shift the noise to align it with the speech signal
alignedNoise = shiftData(noise,-lags(I));

% Subtract shifted noise from speech signal
csSignal = signal - alignedNoise;

    function shiftedData = shiftData(data,lag)
        %Shifts DATA to the right by LAG samples
        if lag==0
            shiftedData=data;
            return
        elseif lag<0
            shiftedData = zeros(size(data));
            shiftedData(1:length(data)+lag) = data((1-lag):length(data));
            shiftedData((length(data)+lag+1):length(data)) = 0;
        else
            shiftedData = zeros(size(data));
            shiftedData(lag+1:length(data)) = data(1:(length(data)-lag));
        end
    end

end

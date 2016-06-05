% Demo.m - Demo for showing noise cancellation using the methods described
% in the paper below.
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

%% Set demoType to "full" (default) or "short"
clear;
tic
demoType = 'full'; % May take up to 10 minutes or more, depending on processor speed. 
% "full" reproduces results from paper: six phrases and one synthesized signal

%demoType = 'short'; % One phrase and one synthesized signal
if(strcmp(demoType,'full'))
    numSignals = 6;
    coeffArray = 0:0.05:1; % Coefficients to try for SNG
else
    numSignals = 1;
    coeffArray = 0:0.05:1;
end

%% Load demo data
addpath subtightplot;
load data/Demo_Full.mat;
signalMat = signalMat(:,1:numSignals);
phrases = phrases(1:numSignals);
[sigLength ~] = size(signalMat);
numCoeffs = length(coeffArray);

%% Add a synthetic signal as ground truth for distortion metric
synSignal = zeros(size(noise3));
synFreqArray = 300:50:3000; % Voice frequency range
x = (1:length(synSignal))';
% Add multiple sine waves from voice frequency range
for freq = synFreqArray
    randPhase = rand(1);
    radFreq = freq/Fs;
    synSignal = synSignal + sin(2*pi*(radFreq*x + randPhase));
end

% Normalize so that power of synthesized signal is equal to power of added
% noise signal (i.e., noise2)
rmsSyn = sqrt(sum(synSignal.^2)/length(synSignal));
rmsNoise2 = sqrt(sum(noise2.^2)/length(noise2));
synSignal = synSignal/rmsSyn*rmsNoise2;
% Remove first 2 seconds of signal (so that there will be a noise only period)
synSignal(1:(2*Fs)) = 0;
% Add noise
synSignalAndNoise = synSignal + noise2;

% Add synthesized signal to signalMat and add to phrases
signalMat = [signalMat synSignalAndNoise];
numSignals = 1+numSignals;
phrases{end+1} = '(Synthesized Signal)';

%% Setup arrays
csSignalArray = cell(numSignals,1); csNoiseArray = cell(numSignals,1); csSupArray = zeros(numSignals,1);
sngSignalArray = cell(numSignals, numCoeffs); sngNoiseArray = cell(numSignals, numCoeffs); sngSupArray = zeros(numSignals, numCoeffs); sngSupBest = zeros(numSignals, 1); sngSignalBest = cell(numSignals, 1); sngCoeffBest = zeros(numSignals, 1);
csSngSignalArray = cell(numSignals, numCoeffs); csSngNoiseArray = cell(numSignals, numCoeffs); csSngSupArray = zeros(numSignals, numCoeffs); csSngSupBest = zeros(numSignals, 1); csSngSignalBest = cell(numSignals, 1); csSngCoeffBest = zeros(numSignals, 1);

for sigIdx = 1:numSignals
    signal = signalMat(:,sigIdx);
    fprintf('\nCanceling noise for phrase "%s"\n', phrases{sigIdx});
    %% Correlation Subtraction (CS)
    csSignalArray{sigIdx} = correlationsubtraction(signal, noise);
    csNoiseArray{sigIdx} = correlationsubtraction(noise2, noise);
    csSupArray(sigIdx) = computenoisereduction('CS', signal, csSignalArray{sigIdx}, noise, csNoiseArray{sigIdx});

    %% Spectral Noise Gating (SNG) with SoX
    % Try multiple noise reduction coefficients
    for i = 1:length(coeffArray)
        sngSignalArray{sigIdx,i} = spectralnoisegating(signal, [], Fs, coeffArray(i));
        sngNoiseArray{sigIdx,i} = spectralnoisegating(noise2, [], Fs, coeffArray(i));
        sngSupArray(sigIdx,i) = computenoisereduction(sprintf('SNG (c = %g)', coeffArray(i)), signal, sngSignalArray{sigIdx,i}, noise, sngNoiseArray{sigIdx,i});
    end
    % Choose signal with the best SNR
    [sngSupBest(sigIdx) bestIdx] = max(sngSupArray(sigIdx,:));
    sngSignalBest{sigIdx} = sngSignalArray{sigIdx, bestIdx};
    sngCoeffBest(sigIdx) = coeffArray(bestIdx);

    %% CS then SNG
    % Try multiple noise reduction coefficients
    for i = 1:length(coeffArray)
        csSngSignalArray{sigIdx,i}= spectralnoisegating(csSignalArray{sigIdx}, [], Fs, coeffArray(i));
        csSngNoiseArray{sigIdx,i} = spectralnoisegating(csNoiseArray{sigIdx}, [], Fs, coeffArray(i));
        csSngSupArray(sigIdx,i) = computenoisereduction(sprintf('CS+SNG (c = %g)', coeffArray(i)), signal, csSngSignalArray{sigIdx,i}, noise, csSngNoiseArray{sigIdx,i});
    end
    % Choose signal with the best SNR
    [csSngSupBest(sigIdx) bestIdx] = max(csSngSupArray(sigIdx,:));
    csSngSignalBest{sigIdx} = csSngSignalArray{sigIdx, bestIdx};
    csSngCoeffBest(sigIdx) = coeffArray(bestIdx);
end

%% Output a few suppression numbers
bestCsSngArray = max(csSngSupArray,[],2);   % Choose suppression with best coefficients
fprintf('\nMinimum CS+SNG suppression for all signals (including synthesized):   %d dB\n',...
        floor(min(bestCsSngArray)));

fprintf('\nCS+SNG suppression for phrase "%s":   %d dB\n',...
        phrases{1}, floor(bestCsSngArray(1)));

fprintf('\nCS+SNG suppression for all other phrases (including synthesized):  %d-%d dB\n', ...
    floor( min(bestCsSngArray(2:end)) ), ceil( max(bestCsSngArray(2:end)) ) );

%% Plot suppression vs. noise reduction coefficient for SNG methods
figure;
plot(coeffArray, [csSngSupArray; sngSupArray], '.-');
hline = findobj(gcf, 'type', 'line');
set(hline((1:numSignals)),'LineStyle',':');
strLegends = strcat([strcat('CS+SNG: "', phrases); strcat('SNG: "', phrases)], '"');
strLegends = reshape(strLegends, length(strLegends(:)), 1); % Fix shape
if(length(strLegends) <= 4)
    legend(strLegends, 'Location', 'Best');
else
    legend(strLegends, 'Location', 'BestOutside');
end

xlabel('Noise Reduction Coefficient'); ylabel('Noise Suppresion (in dB)');
xlim([0 1]); ylim([0 110]);

%% Plot best suppression
figure;
allSup = [csSupArray sngSupBest csSngSupBest]';
labels = {'CS' sprintf('SNG (c_1 = %g)', sngCoeffBest(1)) sprintf('CS+SNG (c_1 = %g)', csSngCoeffBest(1))};
bar(allSup);
legend(phrases, 'Location', 'SE');
set(gca,'XTickLabel',labels'); ylabel('Noise Suppression (in dB)');
set(gca,'YGrid','on');

%% Compute distortion measure on synthesized signal
synIdx = numSignals; % Last signal is synthesized signal
L = length(signal);
NFFT = 2^nextpow2(L); % Next power of 2 from length of y
freqs = Fs/2*linspace(0,1,NFFT/2+1);

signalTitleArray = {'Synthesized', 'Original', 'CS', 'SNG', 'CS+SNG'};
fftSignalArray = cell(5,1);

% Compute FFT of original synthesized signal and synthesized + noise
fftTemp = fft(synSignal,NFFT)/L;
fftSignalArray{1} = 2*abs(fftTemp(1:NFFT/2+1)); % Single sided amplitude
fftTemp = fft(signalMat(:,end),NFFT)/L;
fftSignalArray{2} = 2*abs(fftTemp(1:NFFT/2+1)); % Single sided amplitude

% Compute distortion for all 3 best filtered signals
bestSignals = {csSignalArray{synIdx}, sngSignalBest{synIdx}, csSngSignalBest{synIdx}};
distortion = zeros(3,1);
for i = 1:3
    filtSynSignal = padarray(bestSignals{i}, length(synSignal)-length(bestSignals{i}), 'post');
    
    % Calculate L2 Norm distortion
    distortion(i) = norm(synSignal - filtSynSignal);
    
    % Calculate FFT for next figure
    fftTemp = fft(filtSynSignal,NFFT)/L;
    fftSignalArray{i+2} = 2*abs(fftTemp(1:NFFT/2+1)); % Single sided amplitude
end
bar(distortion);
title('L^2 Norm of Signal Difference (Time domain)');
labels = {'CS' sprintf('SNG (c_syn = %g)', sngCoeffBest(synIdx)) sprintf('CS+SNG (c_syn = %g)', csSngCoeffBest(synIdx))}; 
set(gca,'XTickLabel',labels'); set(gca,'YGrid','on');

%% Plot FFT of synthesized signal
figure;
subplot = @(m,n,p) subtightplot(m, n, p, [0.01 0.01], [0.10 0.10], [0.15 0.04]);
nPlots = length(signalTitleArray);
for i = 1:nPlots
    subplot(nPlots,1,i); 
    subset = 1:round(length(freqs)/7);
    plot(freqs(subset),fftSignalArray{i}(subset));
    if(i == 1); title('Single-Sided Amplitude Spectrum |X(f)|'); end;
    xlabel('Frequency (Hz)');
    set(gca, 'yticklabel', signalTitleArray{i}, 'ytick', [0.018], 'xticklabel', []);
    ylim([0 0.04]);
end
clear subplot;

%% Plot first signal
figure;
subplot = @(m,n,p) subtightplot(m, n, p, [0.01 0.01], [0.10 0.08], [0.12 0.04]);

signalArray = {signalMat(:,1), csSignalArray{1}, sngSignalBest{1}, csSngSignalBest{1}};
titleArray = {'Original', 'CS', 'SNG', 'CS+SNG'};
nPlots = length(titleArray);
for i = 1:nPlots
    subplot(nPlots,1,i); 
    plot((1:length(signalArray{i}))/Fs,signalArray{i});    
    if(i == 1); title(['Phrase "' phrases{1} '"']); end;
    set(gca,'xlim',[0 4]);
    set(gca,'ylim',[-0.8 0.8]);
    set(gca, 'yticklabel', titleArray{i}, 'ytick', [0], 'xticklabel', []);
end
set(gca, 'xticklabelmode', 'auto'); % Reset xticklabel for last
set(gca,'xtick',0:4);
xlabel('Time (s)');
clear subplot;

%% Plot FFT of first signal
figure;
subplot = @(m,n,p) subtightplot(m, n, p, [0.01 0.01], [0.10 0.08], [0.12 0.04]);

L = length(signal);
NFFT = 2^nextpow2(L); % Next power of 2 from length of y
freqs = Fs/2*linspace(0,1,NFFT/2+1);
signalArray = {signalMat(:,1), csSignalArray{1}, sngSignalBest{1}, csSngSignalBest{1}};
titleArray = {'Original', 'CS', 'SNG', 'CS+SNG'};
nPlots = length(titleArray);
for i = 1:nPlots
    curSignal = padarray(signalArray{i}, L-length(signalArray{i}), 'post');
    subplot(nPlots,1,i); 
    X = fft(curSignal,NFFT)/L;
    Xabs = 2*abs(X(1:NFFT/2+1));
    subset = 1:round(length(freqs)/6);
    plot(freqs(subset),Xabs(subset));
    if(i == 1); title(['Single Sided Amplitude Spectrum |X(f)| of Phrase "' phrases{1} '"']); end;
    set(gca, 'yticklabel', titleArray{i}, 'ytick', [0.01], 'xticklabel', []);
    ylim([0 0.02]);
end
xlabel('Frequency (Hz)');
clear subplot;

%% Spectrogram of first signal
figure;
subplot = @(m,n,p) subtightplot(m, n, p, [0.01 0.01], [0.12 0.08], [0.12 0.04]);

signalArray = {signalMat(:,1), csSignalArray{1}, sngSignalBest{1}, csSngSignalBest{1}};
titleArray = {'Original', 'CS', 'SNG', 'CS+SNG'};
nPlots = length(titleArray);
for i = 1:nPlots
    subplot(nPlots,1,i); 
    curSignal = padarray(signalArray{i}, L-length(signalArray{i}), 'post');

    [S,F,T,P]=spectrogram(curSignal,512);
    dbPower = 10*log10(flipud(P));
    dbPower(isinf(dbPower)) = NaN; % Caused by padding vectors to the same length

    imagesc(dbPower);
    if(i == 1); title(['Spectrogram of Phrase "' phrases{1} '"']); end;
    colormap(gray);
    set(gca,'xlim',[0 4/256*Fs]);
    set(gca, 'yticklabel', titleArray{i}, 'ytick', [128], 'xtick', (1:4)/256*Fs, 'xticklabel', []);
    caxis([-180 0]);
    colorbar;
end
set(gca, 'xtick', (1:4)/256*Fs, 'xticklabel', 1:4);
xlabel('Time (s)');

%% Plot final figure of distortion over noise suppression for each method for synthesized signal
figure
ph=plot([csSupArray(synIdx) sngSupBest(synIdx) csSngSupBest(synIdx)],distortion,'ok');
set(ph,'markersize',10,'markerfacecolor','k')
axis([0 120 0 25])
hold on;
% Plot previous results as dotted lines
plot([21 21],[0 25],'k--');
plot([28 28],[0 25],'k--');
xlabel('Noise power suppression (dB)')
ylabel('Distortion*')

%% Play first filtered signal
sound(signalMat(:,1), Fs);
sound(csSignalArray{1}, Fs);
sound(sngSignalBest{1}, Fs);
sound(csSngSignalBest{1}, Fs);
toc

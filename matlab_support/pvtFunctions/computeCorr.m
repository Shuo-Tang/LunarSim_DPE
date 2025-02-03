function [rMax, caf, fdBin] = computeCorr(constellation, signalRx, localCode, settings, freqBin)
% COMPUTECORR Compute correlation result givem the signal and local code
%
% INPUTS:
%   - signalRx: (vector) Received signal (complex-valued).
%   - localCode: (vector) Local code (binary values, e.g., 1 or -1).
%
% OUTPUTS:
%   - rMax: (vector) the maximum correaltion sequence among all frequency bins.
%   - fd: (matrix) 2d correaltion grid based on frequency bins and shift samples .
    
    % With frequency Bins
    % Frequency bins
    % freqBin = settings.frequencyBin; % Doppler frequency bins (Hz)
    numFreqBins = length(freqBin);
    numberOfSamples = length(signalRx);
    if constellation ~= 'E'
        caf = zeros(numFreqBins, settings.samplesPerMs);
    else
        caf = zeros(numFreqBins, settings.samplesPerMs*4);
    end
    %% step computing 
    % % Initialize variables for maximum CAF
    % maxCaf = -1;
    % maxFreqIndex = -1;
    % 
    % % rPlot = zeros(numberOfSamples, numFreqBins);
    % % Loop through each frequency bin
    % for iFreq = 1:numFreqBins
    %     localSignal = zeros(1, numberOfSamples);
    %     % Get the Doppler frequency
    %     fd = freqBin(iFreq);
    %     for iMs = 1:settings.coherentIntegrationTime
    %         startIdx = 1 + (iMs - 1) * settings.samplesPerMs;
    %         endIdx = iMs * settings.samplesPerMs;
    %         ii = 1:settings.samplesPerMs;
    %         % Compute phase values
    %         phase = 2.0 * pi * fd * (ii - 1) / settings.samplingFrequency;
    % 
    %         % Compute the local signal with phase adjustment
    %         localSignal(startIdx: endIdx) = localCode(startIdx: endIdx) .* exp(1i * phase);
    %     end
    % 
    %     % Compute the cross-correlation using FFT-based convolution
    %     r = abs(ifft(fft(signalRx) .* conj(fft(localSignal)))).^2;
    %     if constellation ~= 'E'
    %         caf(iFreq, :) = r(1:settings.samplesPerMs);
    %     else
    %         caf(iFreq, :) = r(1:settings.samplesPerMs*4);
    %     end
    % 
    %     % find maximum correlation sequence
    %     maxR = max(r);
    %     if maxR > maxCaf
    %         maxFreqIndex = iFreq;
    %         maxCaf = maxR;
    %     end
    % end
    % rMax = caf(maxFreqIndex,:);

    %% parallel computing
    % Initialize variables for maximum CAF
    rMaxFreq = zeros(length(freqBin), 1);
    % rPlot = zeros(numberOfSamples, numFreqBins);
    % Loop through each frequency bin
    parfor iFreq = 1:numFreqBins
        localSignal = zeros(1, numberOfSamples);
        % Get the Doppler frequency
        fd = freqBin(iFreq);
        for iMs = 1:settings.coherentIntegrationTime
            startIdx = 1 + (iMs - 1) * settings.samplesPerMs;
            endIdx = iMs * settings.samplesPerMs;
            ii = 1:settings.samplesPerMs;
            % Compute phase values
            phase = 2.0 * pi * fd * (ii - 1) / settings.samplingFrequency;

            % Compute the local signal with phase adjustment
            localSignal(startIdx: endIdx) = localCode(startIdx: endIdx) .* exp(1i * phase);
        end

        % Compute the cross-correlation using FFT-based convolution
        r = abs(ifft(fft(signalRx) .* conj(fft(localSignal)))).^2;
        rMaxFreq(iFreq) = max(r);
        if constellation ~= 'E'
            caf(iFreq, :) = r(1:settings.samplesPerMs);
        else
            caf(iFreq, :) = r(1:settings.samplesPerMs*4);
        end   
    end
    % find maximum correlation sequence
    [~, maxFreqIndex] = max(rMaxFreq);

    rMax = caf(maxFreqIndex,:);
    %% set next frequency bin

    fdBin = freqBin(maxFreqIndex)-3*500:500:freqBin(maxFreqIndex)+4*500;

    % % Visualize CAF
    % cafPlot2SP(rPlot, freqBin, 1:numberOfSamples);
    % % Extract the Doppler frequency and pseudorange
    % fd = freqBin(maxFreqIndex);
    % pr = (maxSampleIndex - 1) / settings.samplingFrequency * constants.c; % Convert sample index to time, then to meters

end
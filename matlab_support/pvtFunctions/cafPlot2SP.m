function cafPlot2SP(r, Doppler, sampleIndices)
% PLOT3DCAF Plots the cross-ambiguity function (CAF) in 3D.
%
% INPUTS:
%   - r: (matrix) CAF values (numberOfSamples x numFreqBins).
%   - Doppler: (vector) Doppler frequency bins (Hz).
%   - sampleIndices: (vector) Number of samples bins.
%
% The plot will display the CAF with Doppler frequencies on one axis,
% sample indices on another axis, and the CAF magnitude as the vertical axis.

    % Create meshgrid for Doppler and sample indices
    [dopplerGrid, sampleGrid] = meshgrid(Doppler, sampleIndices);

    % Plot the 3D surface
    figure;
    surf(dopplerGrid, sampleGrid, abs(r));
    % shading interp; % Smooth shading for better visualization
    colorbar; % Add a color bar to indicate magnitude
    xlabel('Doppler Frequency (Hz)');
    ylabel('Sample Index');
    zlabel('CAF Magnitude');
    title('Cross-Ambiguity Function (CAF)');
    view(0, 0); % Set a good viewing angle
end

function progressBar(iUpdate, numOfUpdates)
    % Function to display a single-line progress bar in the command window.
    % Input:
    %   iUpdate       - Current iteration (1-based index)
    %   numOfUpdates  - Total number of iterations

    % Total width of the progress bar
    totalBarWidth = 50;

    % Calculate progress percentage
    progress = iUpdate / numOfUpdates * 100;

    % Calculate the number of '#' and spaces
    numHashes = floor(progress / 100 * totalBarWidth);
    numSpaces = totalBarWidth - numHashes;

    % Construct the progress bar string
    progressBar = sprintf('Progress: [%s%s] %3.0f%%', ...
        repmat('#', 1, numHashes), repmat(' ', 1, numSpaces), progress);

    % Use '\r' to overwrite the previous line in the command window
    fprintf('\r%s', progressBar);

end
function EstRxPVT = conv2steps(r,UserPosition,SatPosition, settings, constants)

%% load configuration file
fs = settings.samplingFrequency;
c = constants.c;
num2stepsIterations = settings.rlsIter;
numSV = size(r, 1);
%% memory allocation
EstRange=zeros(1,numSV);


RefPos=UserPosition;
EstRxPVT=RefPos;


%Estimate time delays
[~, maxPos] = max(r,[],2);
maxPos=maxPos-1;



EstFracDelay=maxPos/fs;
EstFracRange=EstFracDelay * c;

% Loop over iterations.
for kIterations                                     =   1 : num2stepsIterations
    
    
    % Loop over satellites.
    for kSV                                         =   1 : numSV
        EstRange(kSV)                               =   norm(SatPosition(kSV,:) - EstRxPVT(1:3));
        
        numH                                        =   SatPosition(kSV, :) - EstRxPVT(1:3);
        denH                                        =   norm(numH);
        H(kSV, 1:3)                      =   - numH / denH;
        %          H(kSV,4)=1;
        
    end
    
    % corrP                                                   =   (EstFracRange - EstRange') / c;
    % corrP_noAmbg                                            =   wrap(rem(corrP, 1e-3), 0.5e-3);
    % corrFracPseudorange                                     =   corrP_noAmbg * c;
    fracRange = rem(EstRange / c, 1e-3) * c;
    corrP     = EstFracRange - fracRange';

    % deltaPVT                                                =   ((H' * H) \ H') * corrFracPseudorange;
    deltaPVT                                                =   ((H' * H) \ H') * corrP;
    EstRxPVT                                         =   EstRxPVT + deltaPVT.';
end

PosErrLS=norm(EstRxPVT(1:3)-UserPosition);
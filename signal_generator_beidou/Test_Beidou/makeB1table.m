function B1CodesTable = makeB1table(settings)
%Function generates B1 codes for all 37 satellites based on the settings
%provided in the structure "settings". The codes are digitized at the
%sampling frequency specified in the settings structure.
%One row in the "B1CodesTable" is one B1 code. The row number is the PRN
%number of the B1 code.
%
%caCodesTable = makeCaTable(settings)
%
%   Inputs:
%       settings        - receiver settings
%   Outputs:
%       B1CodesTable    - an array of arrays (matrix) containing B1 codes
%                       for all satellite PRN-s

% Find number of samples per spreading code ----------------------------
samplesPerCode = round(settings.samplingFreq / (settings.codeFreqBasis / settings.codeLength));

% Prepare the output matrix to speed up function -----------------------
B1CodesTable = zeros(37, samplesPerCode);

  % Find time constants --------------------------------------------------
ts = 1/settings.samplingFreq;   % Sampling period in sec
tc = 1/settings.codeFreqBasis;  % B1 chip period in sec

  % For all satellite PRN-s ...
for No = 1:37
    % Generate B1 code for given PRN -----------------------------------
    B1Code = generateB1code(No);
 
    % Digitizing =======================================================
    
    % Make index array to read B1 code values -------------------------
    % The length of the index array depends on the sampling frequency -
    % number of samples per millisecond (because one B1 code period is one
    % millisecond).
    codeValueIndex = ceil((ts * (1:samplesPerCode)) / tc);
    
    % Correct the last index (due to number rounding issues) -----------
    codeValueIndex(end) = 2046;
    % Make the digitized version of the B1 code -----------------------
    % The "upsampled" code is made by selecting values form the B1 code
    % chip array (caCode) for the time instances of each sample.
    B1CodesTable(No, :) = B1Code(codeValueIndex);
end, % for PRN = 1:37
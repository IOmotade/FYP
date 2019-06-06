% Omotade Iluromi, GROUP (EE4), 2019, Imperial College.
% 26/05/2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate readable format for measured values along with units and
% prefixes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs
% x (RxC Double) = Values
% unit (String) = Units of Measurement
% nd (Integer) = Number of decimal places for output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Outputs
% z (RxC String) = Values of x with appropriate S.I. prefixes and defined
% unit
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function z = fUnits(x, unit, nd)
if ~exist('unit', 'var')
    unit = '';
end
if ~exist('nd', 'var')
    nd = 2;
end
nd = repmat(nd, size(x));
unit = string(unit);
sameSizeCheck = (numel(x) == numel(unit));
if ~sameSizeCheck
    unit = repmat(unit, size(x));
end
z = arrayfun(@fUnit, x, unit, nd);
end

function z = fUnit(xPwr, unit, nd)
prefix.p = ['y', 'z', 'a', 'f', 'p', 'n', 'u', 'm', ' ',...
    'k', 'M', 'G', 'T', 'P', 'E', 'Z', 'Y']; %Must be single characters
prefix.pwr = (-24:3:24);
%% Extract Signs
[sign, xPwr] = deal(xPwr/abs(xPwr), abs(xPwr));
sign = char((' ')*(sign==1) + ('-')*(sign==-1));
if xPwr==0
    z = "0";
    return
end
%% Extract Index Powers
idxPwr = log10(xPwr);
xPwr = idxPwr - round(idxPwr);
idxPwr = round(idxPwr);

%% Convert to SI units
siPwr = 3*floor(idxPwr/3);
xPwr = xPwr+idxPwr-siPwr;
prefixIdx = (siPwr == prefix.pwr);
prfx = prefix.p(prefixIdx);

%% Check if it could be found in the units
if(isempty(prfx))
    xPwr = xPwr+(siPwr - prefix.pwr(end));
    prfx = prefix.p(end);
end

%% Present in SI Units
x = 10^xPwr;
txt = sprintf('%s%d%s', '%s%2.', nd, 'f%s%s');
z = string(sprintf(txt, sign, x, prfx, unit));
end
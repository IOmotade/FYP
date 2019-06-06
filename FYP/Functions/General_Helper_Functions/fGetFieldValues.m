% Omotade Iluromi, GROUP (EE4), 2019, Imperial College.
% 26/05/2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Converts a simulated circuit from the time domain to frequency domain
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs
% structure (Struct) = Structure object with desired field
% varargin (1xN String) = List of fields till desired fields e.g.
% X.a.b.c => fx(X, 'a', 'b', 'c')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Outputs
% val (AxBxCxD Type) = Matrix containing entries of desired field in
% structure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function val = fGetFieldValues(structure, varargin)

[rowSize, colSize] = size(structure);

scheme = length(varargin);
for rowIdx = 1:rowSize
    for colIdx = 1:colSize
        val(rowIdx, colIdx, :, :) = fGetField(structure, rowIdx, colIdx, varargin, scheme);
    end
end
val = squeeze(val);
end

function f = fGetField(structure, rowIdx, colIdx, var, scheme)
if scheme ==1
    f = getfield(structure, {rowIdx, colIdx}, char(var{1}));
    return
end

if scheme ==2
    f = getfield(structure, {rowIdx, colIdx}, char(var{1}), char(var{2}));
    return
end

if scheme ==3
    f = getfield(structure, {rowIdx, colIdx}, char(var{1}), char(var{2}), char(var{3}));
    return
end

if scheme ==4
    f = getfield(structure, {rowIdx, colIdx}, char(var{1}), char(var{2}), char(var{3}), char(var{4}));
    return
end

end
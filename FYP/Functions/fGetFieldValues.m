function val = fGetFieldValues(structure, varargin)

% field = char(field);
[rowSize, colSize] = size(structure);

scheme = length(varargin);
for rowIdx = 1:rowSize
    for colIdx = 1:colSize
        val(rowIdx, colIdx, :, :) = fGetField(structure, rowIdx, colIdx, varargin, scheme);
        %         val(rowIdx, colIdx, :, :) = getfield(structure, {rowIdx, colIdx}, field, 'value');
        %         val(rowIdx, colIdx, :, :) = getfield(structure, {rowIdx, colIdx}, field);
    end
end
val = squeeze(val);
end

function f = fGetField(structure, rowIdx, colIdx, var, scheme)
if scheme ==1
    f = getfield(structure, {rowIdx, colIdx}, char(var{1}));
end

if scheme ==2
    f = getfield(structure, {rowIdx, colIdx}, char(var{1}), char(var{2}));
end

if scheme ==3
    f = getfield(structure, {rowIdx, colIdx}, char(var{1}), char(var{2}), char(var{3}));
end

if scheme ==4
    f = getfield(structure, {rowIdx, colIdx}, char(var{1}), char(var{2}), char(var{3}), char(var{4}));
    
end

end
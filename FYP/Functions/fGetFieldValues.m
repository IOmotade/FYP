function val = fGetFieldValues(structure, field)

field = char(field);
[rowSize, colSize] = size(structure);

for rowIdx = 1:rowSize
    for colIdx = 1:colSize
        val(rowIdx, colIdx, :, :) = getfield(structure, {rowIdx, colIdx}, field, 'value');
    end
end
val = squeeze(val);
end
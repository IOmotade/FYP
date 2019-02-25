function fSetFieldValues(structure, field, val)

field = char(field);
% [tSize, rowSize, colSize] = size(val);
[tSize, ~, ~] = size(val);

% for rowIdx = 1:rowSize
%     for colIdx = 1:colSize
for tIdx = 1:tSize
    structure = setfield(structure, {tIdx}, field, 'value', val(tIdx, :, :));
end
% end

end
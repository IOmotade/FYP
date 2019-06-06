function fieldExist = fFieldExist (inStruct, fieldName)
% inStruct is the name of the structure or an array of structures to search
% fieldName is the name of the field for which the function searches
fieldExist = 0;
f = fieldnames(inStruct(1));
for i=1:length(f)
    if(strcmp(f{i},strtrim(fieldName)))
        fieldExist = 1;
        return;
    elseif isstruct(inStruct(1).(f{i}))
        fieldExist = fFieldExist(inStruct(1).(f{i}), fieldName);
        if fieldExist
            return;
        end
    end
end
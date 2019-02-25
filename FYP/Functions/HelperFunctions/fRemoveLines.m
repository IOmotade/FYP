function fRemoveLines(fileLoc, lineChars)
if size(fileLoc, 1) ~= size(lineChars, 1)
    fileLoc = repmat(fileLoc, size(lineChars));
end
arrayfun(@fRemoveLine, fileLoc, lineChars);
end

function fRemoveLine(fileLoc, lineChars)

[filePath, name, ext] = fileparts(fileLoc);
tmpFileLoc = strcat(filePath, '/tmp', name, ext);
expr = char(strcat(lineChars, '(\w*)'));

%Make temporary copy of file
copyfile(char(fileLoc), char(tmpFileLoc))

%Open Files
fid_in  = fopen(tmpFileLoc);
fid_out = fopen(fileLoc,'w+');

text_in = fgetl(fid_in);
while ischar(text_in)
    %Check if an empty line
    valid = ~isempty(text_in);
    
    %Remove unneccessary lines
    strtIdx = regexp(text_in,expr, 'ONCE');
    valid = valid & isempty(strtIdx);
    if  valid
        fprintf(fid_out,'%s\n',text_in);
    end
    text_in = fgetl(fid_in);
end

fclose(fid_in);
fclose(fid_out);

delete(tmpFileLoc)

end
function fRemoveExtraWhitespace(fileLoc)

[filePath, name, ext] = fileparts(fileLoc);
tmpFileLoc = strcat(filePath, '/tmp', name, ext);

%Make temporary copy of file
copyfile(char(fileLoc), char(tmpFileLoc))

%Open Files
fid_in  = fopen(tmpFileLoc);
fid_out = fopen(fileLoc,'w+');

text_in = fgetl(fid_in);
while ischar(text_in)
    % Clear leading whitespace
    text_out = strtrim(text_in);
    
    % Remove other whitespaces
    text_out = regexprep(text_out,'(\t|\s+)',', ');
    fprintf(fid_out,'%s\n',text_out);
    text_in = fgetl(fid_in);
end

fclose(fid_in);
fclose(fid_out);

delete(tmpFileLoc)

end
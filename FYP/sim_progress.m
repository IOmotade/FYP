%% Get current sims progress
rst
% copyfile('sim.log', 'tmp_log.txt');
unix('cat sim.log | col -b > tmp_log.txt');

filename = 'tmp_log.txt';
fid = fopen(filename, 'r');
[lineIdx, numLinesToSkip] = deal(0, 9);
x = '';


tline = fgetl(fid);
while ischar(tline) && lineIdx<numLinesToSkip
    lineIdx = lineIdx + 1;
    %     x = sprintf('%s, %s', x, tline);
    tline = fgetl(fid);
end

while ischar(tline)
    %     tline(tline==' ') = [];
    tline(tline==',') = [];
    x = sprintf("%s\n%s", x, tline);
    tline = fgetl(fid);
end

x = splitlines(x);

id = " ";
id_lngth = 12;
for idx= 1:length(x)
    tmp = char(x(idx));
    if length(tmp)>id_lngth
        tmp_id = tmp(1:id_lngth);
        if any(id==string(tmp_id))
            txt(id==string(tmp_id)) = x(idx);
        else
            id = [id; string(tmp_id)];
        end
    end
end

% [id, order] = sort(id);
% txt = txt(order);
txt(isempty(txt))=[];
for idx = 1:length(txt)
    if ~ismissing(txt(idx))
        disp(txt(idx))
    end
end


fclose(fid);

delete(filename)
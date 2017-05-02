function [metadata] = read(fname, traces)

fileID = fopen(fname, 'rb');

metadata.traces = fread(fileID, 1,'uint32');
metadata.key = fread(fileID, 16, 'uint8');
for i = 1:traces
    metadata.P{i} = fread(fileID, 16, 'uint8');
    metadata.C{i} = fread(fileID, 16, 'uint8');
end
fclose(fileID);
end
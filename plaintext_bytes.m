function [result] = plaintext_bytes(struct, traces)
result = [];
for i = 1:traces
    s = struct.P{i};
%     a = s{1}(byte);
    result = [result; s];
end
end
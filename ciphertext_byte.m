function [result] = ciphertext_byte(struct, byte, traces)
result = [];
for i = 1:traces
    s = struct.C{i}(byte);
%     a = s{1}(byte);
    result = [result; s];
end
end
function [A] = load_data(fname, nr)
A = []; 
for i=1:nr
target_file = strcat(fname, int2str(i));
load(target_file);
% A = data;
A = [A; data];
end
A = double(A);
end
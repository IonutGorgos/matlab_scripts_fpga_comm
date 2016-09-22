function H = hamming_weight(X)
% HAMMING_WEIGHT Computes the Hamming weight
%   [H] = hamming_weight(X)
%   returns the Hamming weight of the elements in X.

%% Initialize and check parameters
[m, n] = size(X);
H = zeros(m,n);

%% Compute the Hamming weight of each element in X
for i = 1:m
    for j=1:n
        valbin = dec2bin(X(i,j));
        for k=1:length(valbin)
            H(i,j) = H(i,j) + (valbin(k) == '1');
        end
    end
end

end
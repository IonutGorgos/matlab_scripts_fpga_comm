function [comb] = combine_aes_sbox(plaintext, key, cf_params)
%COMBINE_AES_SBOX Returns the AES S-box of the XOR of plaintext and key
%   [comb] = COMBINE_AES_SBOX(plaintext, key, s_box)
%   returns the AES S-box output of the XOR between a plaintext and key
%   byte (or bytes).
%
%   plaintext and key should be vectors of length nr_values to be combined
%   through XOR and then applied to the AES S-box.
%
%   cf_params is a structure having optional parameters for this method.
%   For this particular method, the only useful parameter is:
%
%   - 's_box': which provides a method or vector that
%   computes/returns the AES S-box. If this is not provided, or empty ([]),
%   then a default version is used. Providing it might improve speed. If
%   provided, the input values should be considered as starting from "1",
%   to provide compatibility with MATLAB vectors. Note that this might
%   introduce some errors if dealing with "uint8" values, so check your
%   use/implementation with this code.

%% Check and initialize data
nr_values = length(plaintext);
if nr_values ~= length(key)
    if length(key) == 1
        key = repmat(key, nr_values, 1);
    else
        error('plaintext and key have different lengths');
    end
end
if nargin < 3 || isempty(cf_params) || ~isfield(cf_params, 's_box')
    [s_box, ~] = s_box_gen;
else
    s_box = cf_params.s_box;
end

bx = double(bitxor(plaintext, key));
comb = s_box(bx+1);
                              
end
    

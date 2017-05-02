%AES_DEMO  Demonstration of AES-components.
%
%   AES_DEMO
%   runs a demonstration of all components of 
%   the Advanced Encryption Standard (AES) toolbox.
%
%   In the initialization step the S-boxes, the round constants,
%   and the polynomial matrices are created and
%   an example cipher key is expanded into 
%   the round key schedule.
%   Step two and three finally convert 
%   an example plaintext to ciphertext and back to plaintext.

%   Copyright 2001-2005, J. J. Buchholz, Hochschule Bremen, buchholz@hs-bremen.de

%   Version 1.0     30.05.2001

% Initialization
[s_box, inv_s_box, w, poly_mat, inv_poly_mat] = aes_init;

% Define an arbitrary series of 16 plaintext bytes 
% in hexadecimal (string) representation
% The following two specific plaintexts are used as examples 
% in the AES-Specification (draft)
% plaintext_hex = {'00' '11' '22' '33' '44' '55' '66' '77' ...
%                  '88' '99' 'aa' 'bb' 'cc' 'dd' 'ee' 'ff'};
plaintext_hex = {'3A' 'D4' 'f7' '6d' '66' '2a' '16' '02' ...
                'ed' '33' 'eb' '22' '93' 'd0' 'd8' 'dd'};
% plaintext_hex = {'00' '01' '02' '03' '04' '05' '06' '07' ...
%                  '08' '09' '0a' '0b' '0c' '0d' '0e' '0f'};
% Convert plaintext from hexadecimal (string) to decimal representation
plaintext = hex2dec (plaintext_hex);

% This is the real McCoy.
% Convert the plaintext to ciphertext,
% using the expanded key, the S-box, and the polynomial transformation matrix
% ciphertext = cipher (plaintext, w, s_box, poly_mat, 1);

% Convert the ciphertext back to plaintext
% using the expanded key, the inverse S-box, 
% and the inverse polynomial transformation matrix
re_plaintext = inv_cipher (ciphertext, w, inv_s_box, inv_poly_mat, 1);
function[plaintext,ciphertext,traces,key] = gather_data(fname, fname2,fout, nr, nr_traces) 
struct1 = read(fname, nr_traces);
key = struct1.key;
plaintext = []; 
ciphertext = [];
for i = 1:16
plaintext = [plaintext, plaintext_byte(struct1,i,nr_traces)];
ciphertext = [ciphertext, ciphertext_byte(struct1,i,nr_traces)];
end
traces = load_data(fname2, nr);
save(fout, 'plaintext','ciphertext', 'traces','key', '-v7.3');
end
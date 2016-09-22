function[corr_power] = do_cpa_attack(fname, fname2, byte, nr, traces) 

model.values = 0:255;
model.lmodel = hamming_weight(model.values);
% save_625(3)
struct1 = read(fname, traces);
plaintext_power_on = plaintext_byte(struct1,byte,traces);
data_power_on = load_data(fname2, nr);
data_power_on = adc2mv(data_power_on,50*10^-3,32512);
corr_power = compute_scores_cpa_sbox(data_power_on,plaintext_power_on,'classic',model);
end
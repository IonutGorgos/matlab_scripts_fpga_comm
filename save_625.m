function[] = save_625(nr_file) 
for i=1:nr_file
target_file = strcat('trace_', int2str(i))
load(target_file)
data = data(:,1:625);
Time = Time(1:625);
numSamples = 625;
filename = strcat('data\data_', int2str(i))
save(filename,'data','Time','numSamples')
end
function [xdata] = compute_features_data(data, ...
                                    func_prepare, ...
                                    pp1, pp2, pp3, pp4, pp5)
%COMPUTE_FEATURES Computes the features of a set of power traces.
%   [xdata] = COMPUTE_FEATURES_DATA(data, ...
%                              func_prepare, ...
%                              pp1, pp2, pp3, pp4, pp5)
%   computes some features of a set of given power traces. 
%
%   data should be a matrix of size nr_trials x nr_samples x nr_groups,
%   where:
%   - nr_trials is the number of traces each having nr_samples leakage
%     samples.
%   - nr_groups is the number of groups containing leakage traces
%   This matrix should contain the data to be processed by this method.
%
%   The func_prepare method should transform the data (e.g. compress it)
%   for the further processing. The output data can then be used to test
%   leakage, classification, etc. An exmaple of func_prepare method is
%   prepare_data_pca which uses some PCA parameters to retain only some
%   principal directions of the data. The prototype of func_prepare is:
%   data_out = func_prepare(data_in, pp1, pp2, pp3, ...) where data_in is
%   the raw input trace and pp1, pp2, pp3 are parameters needed by the
%   specific function.
%
%   The output xdata is a matrix of size nr_trials x nr_features x nr_groups,
%   where nr_trials is the number of traces per group as given in data,
%   nr_features is the number of features remaining per trace after
%   applying the func_prepare method and nr_groups is the number of groups
%   specified by data.

%% Check and initialize data
nr_groups = size(data, 3);

%% Process
for k=1:nr_groups
    data_out = func_prepare(data(:,:,k), pp1, pp2, pp3, pp4, pp5);    
    if k==1        
        [m,n] = size(data_out);
        xdata = zeros(m, n, nr_groups);
    end    
    xdata(:,:,k) = data_out;
end
                              
end
    

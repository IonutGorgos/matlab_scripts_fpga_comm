function [sr] = compute_success_rate(results, nr_traces_vec)
%COMPUTE_SUCCESS_RATE Computes the success rates
%   [sr] = COMPUTE_SUCCESS_RATE(results, nr_traces_vec)
%   Computes the success rate for the givn results.
%
%   results should be a cell of results, such as those
%   returned by run_cpa_sbox_scores, for which this method
%   will compute the success rate. This should be a column
%   vector of cells of length nr_iterations, where
%   nr_iterations is the number of iterations that the
%   experiment was carried out.
%
%   nr_traces_vec should be a vector which specifies the
%   number of attack traces used for each experiment.

%% Initialise variables
nr_iter = size(results, 1);
nr_vec = length(nr_traces_vec);
sr = zeros(nr_vec, 1);

%% Compute success rate
for i=1:nr_iter
    for k=1:nr_vec
        rr = results{i}.disc_info.depth.(['group' num2str(k)]);
        if(rr(1) == 1)
            sr(k) = sr(k) + 1;
        end
    end
end
sr = sr ./ nr_iter;

end




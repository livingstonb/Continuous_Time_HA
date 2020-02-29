
clear
load('/home/brian/Documents/temp/income_process.mat')

import HACTLib.random.simulate_markov

markov_options.n_sim = 2e4;
markov_options.t_sim = 1e2;
markov_options.stepsize = 1e-3;
markov_options.normalize = false;

for k = 1:55
    fprintf('Simulating k = %d\n', k)
    markov_options.initial_index = k;

    T = simulate_markov(ytrans, values, markov_options);
    tmp = table2array(T(2,2));
    expectation(k) = tmp{1};
end

expectation = expectation(:);
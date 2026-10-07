function stats = paired_ci(delta)
delta = delta(:);
n = numel(delta);
if n < 2
    error('At least two paired observations are required.');
end
tcrit = 1.98421695150868;
mu = mean(delta);
half = tcrit*std(delta,0)/sqrt(n);
stats = [mu, mu-half, mu+half];
end

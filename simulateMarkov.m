function [ sample, sampleIx ] = simulateMarkov( T, X, Pi )
n = size(Pi, 1);

sampleIx = ones([ T, 1 ]);
sampleIx(1) = floor(n/2);
cumPi = cumsum(Pi, 2);
draws = rand(T, 1);

for tt = 2:T
    for stateToday = 1:n
        if draws(tt) <= cumPi(sampleIx(tt-1), stateToday)
            sampleIx(tt) = stateToday;
            break;
        end
    end
end

sample = X(sampleIx);
end


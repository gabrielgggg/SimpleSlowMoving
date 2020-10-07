function [ X, Pi ] = makeMC(yBar, rho, sgma, n, trimTails)
%
% Construct Markov Chain with states X and transition probabilities Pi
% using AR(1) with mean yBar, autocorrelation rho, and s.d. of innovation
% sgma, with n points of support.
%
% trimTrails = true => discretize truncated Normal
%            = false => assign mass in tails to end points of grid
%

% X = equally spaced over an interval
% [ yBar - 3 s.d, yBar + 3 s.d. ]
X = linspace( yBar - 3.0 * sgma / sqrt(1 - rho^2), ...
    yBar + 3.0 * sgma / sqrt(1 - rho^2), n);

% Pi - transition matrix
Pi = zeros([ n, n ]);
for stYesterday = 1:n
    % work at state stYesterday
    condE = (1-rho) * yBar + rho * X(stYesterday);
    % Interior states...
    for stToday = 2:n-1
        Pi(stYesterday, stToday) = ...
            normcdf((X(stToday) + X(stToday+1))/2, condE, sgma) ...
            - normcdf((X(stToday-1) + X(stToday))/2, condE, sgma);
    end

    if trimTails
        rightOfFirst = (X(1) + X(2))/2;
        shiftOfFirst = rightOfFirst - X(1);
        Pi(stYesterday, 1) = ...
            normcdf(rightOfFirst, condE, sgma) ...
            - normcdf(X(1) - shiftOfFirst, condE, sgma);
        
        leftOfLast = (X(end-1) + X(end))/2;
        shiftOfLast = X(end) - leftOfLast;
        Pi(stYesterday, end) = ...
            normcdf(X(end) + shiftOfLast, condE, sgma) ...
            - normcdf(leftOfLast, condE, sgma);
        
        Pi(stYesterday, :) = Pi(stYesterday, :) ./ sum(Pi(stYesterday, :));
    else
        Pi(stYesterday, 1) = normcdf((X(1) + X(2))/2, condE, sgma);
        Pi(stYesterday, n) = 1 - normcdf((X(n-1) + X(n))/2, condE, sgma);
    end
end

end

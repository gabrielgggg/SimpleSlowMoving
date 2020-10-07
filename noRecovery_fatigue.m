
%
%
% tau0 + tau1( b - bBar) + q (b' - (1-delta) b) = k b
%
% Steady state:
% tau0 + delta q bBar = k bBar
% tau0 = k bBar - delta q bBar
% tau0 = r bBar  if k = delta + r and q = 1
%

clear;
close all;

r = 0.04;
delta = 1.0 / 20.0;
kappa = r + delta;
bBar = 1.0;
tau0 = r * bBar;
tau1 = 0.25;
maxPS = 0.125;

bSz = 501;
zSz = 51;

% b = linspace(-0.25, 0.75, bSz);
b = linspace(-0.5, 2.5, bSz);
[z, zPi] = makeMC(0.0, 0.95, 0.0075, zSz, true);

bInterval = b(2) - b(1);

bPr = zeros(zSz, bSz);
q = zeros(zSz, bSz);
d = zeros(zSz, bSz);

err = 1.0;
iter = 1;
errTol = 1.0e-5;
maxIter = 400;

while err > errTol && iter <= maxIter
  for zIx = 1:zSz
    for bIx = 1:bSz
      bHere = b(bIx);
      laff = q(zIx, :) .* (b - (1.0 - delta) * bHere); %#ok<*PFBNS>
      thold = kappa * bHere - min(maxPS, tau0 + tau1 * (bHere - bBar)) - z(zIx);
      if max(laff) < thold
        if bHere > 0.0
          d(zIx, bIx) = 1;
        else
          d(zIx, bIx) = 0;
        end
        bPr(zIx, bIx) = 0.0;
      else
        d(zIx, bIx) = 0;
        if min(laff - thold) > 0
          bPr(zIx, bIx) = b(1);
        else
          for bPrIx = 1:bSz-1
            if laff(bPrIx) < thold && laff(bPrIx+1) >= thold
              slope = (laff(bPrIx+1) - laff(bPrIx)) / bInterval;
              intercept = laff(bPrIx) - thold - slope * b(bPrIx);
              bPr(zIx, bIx) = -intercept / slope;
              break;
            end
          end
        end
      end
    end % bIx
  end % zIx
  
  newq = zeros(size(q));
  for bIx = 1:bSz
    if b(bIx) <= 0.0
      newq(:, bIx) = 1.0;
    else
      for zIx = 1:zSz
        expTerms = zeros([zIx, 1]);
        for zPrIx = 1:zSz
          bPrPr = bPr(zPrIx, bIx);
          bSegments = floor((bPrPr - b(1)) / bInterval)+1;
          slope = (q(zPrIx, bSegments+1) - q(zPrIx, bSegments)) / bInterval;
          intercept = q(zPrIx, bSegments) - slope * b(bSegments);
          contQ = intercept + slope * bPrPr;
          expTerms(zPrIx) = (1.0 - d(zPrIx, bIx)) * ( kappa + (1.0-delta) * contQ );
        end
        
        newq(zIx, bIx) = dot(zPi(zIx, :), expTerms) / (1.0 + r);
      end
    end % zIx
  end % bIx
  
  err = max(abs(q(:)  - newq(:) ));
  q = newq;
  fprintf('Iter %d with error %f \n', iter, err);
  iter = iter + 1;
end

figure;
subplot(1, 3, 1); plot(b, bPr(1:5:end, :), '-x'); hold on; plot(b, b, 'k--', 'LineWidth', 2); title('b''');
subplot(1, 3, 2); plot(b, d(1:5:end, :)); title('d');
qb = q .* repmat(b, [zSz, 1]);
subplot(1, 3, 3); plot(b, qb(1:5:end, :)); title('q b');

T = 50000;
[ simZ, simZix ] = simulateMarkov( T, z, zPi );
simB = bBar * ones([T, 1]);
simD = zeros([T, 1]);
simQ = zeros([T, 1]);

for tIx = 1:T-1
  bHere = simB(tIx);
  zHere = simZ(tIx);
  laff = q(simZix(tIx), :) .* (b - (1.0 - delta) * bHere); %#ok<*PFBNS>
  thold = kappa * bHere - min(maxPS, tau0 + tau1 * (bHere - bBar)) - zHere;
  if max(laff) < thold
    if bHere > 0.0
      simD(tIx) = 1;
    else
      simD(tIx) = 0;
    end
    simB(tIx+1) = 0.0;
  else
    simD(tIx) = 0;
    if min(laff - thold) > 0
      simB(tIx+1) = b(1);
    else
      for bPrIx = 1:bSz-1
        if laff(bPrIx) < thold && laff(bPrIx+1) >= thold
          slope = (laff(bPrIx+1) - laff(bPrIx)) / bInterval;
          intercept = laff(bPrIx) - thold - slope * b(bPrIx);
          simB(tIx+1) = -intercept / slope;
          break;
        end
      end
    end
  end
end

figure;
subplot(2, 2, 1); histogram(simZ);
subplot(2, 2, 2); histogram(simB);
subplot(2, 2, 3); histogram(simD, 'Normalization', 'probability');

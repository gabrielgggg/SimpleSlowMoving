clear;
close all;

r = 0.04;
delta = 1.0 / 6.0;
kappa = r + delta;
bBar = 0.35;

atB = [ 0.0, bBar, 1.5];
primSurp = [ -0.05, r*bBar, 0.0785 ];
fatigue = @(coefs, bb) coefs(1) +  coefs(2) ./ (1.0 + exp(-coefs(3) * bb));
ccf = fsolve( @(cc) primSurp - fatigue(cc, atB), [-0.1816 0.2631 3.4111]);
prim = @(bb) fatigue(ccf, bb);

figure;
scatter(atB, primSurp);
hold on;
fplot( @(bb) fatigue(ccf, bb), [-0.25, 1.75]);
pause(0.5);

nu = 0.5;
maxRecov = 0.6;
qMin = 0.3;

bSz = 501;
zSz = 51;

% b = linspace(-0.25, 0.75, bSz);
b = linspace(-0.5, 2.0, bSz);
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
      if iter > 30
        laff(q(zIx, :) < qMin) = -1000; % "underwriting standards"
      end
      thold = kappa * bHere - prim(bHere) - z(zIx);
      if max(laff) < thold
        if bHere > 0.0
          d(zIx, bIx) = 1;
          bPr(zIx, bIx) = min(maxRecov, nu * bHere);
        else
          d(zIx, bIx) = 0;
          bPr(zIx, bIx) = 0.0;
        end
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
          
          recovB = min(maxRecov, nu * b(bIx));
          bSegments = floor((recovB - b(1)) / bInterval)+1;
          slope = (q(zPrIx, bSegments+1) - q(zPrIx, bSegments)) / bInterval;
          intercept = q(zPrIx, bSegments) - slope * b(bSegments);
          recovQ = intercept + slope * recovB;
          
          expTerms(zPrIx) = (1.0 - d(zPrIx, bIx)) * ( kappa + (1.0-delta) * contQ ) ...
            + d(zPrIx, bIx) * min(maxRecov / b(bIx), nu) * recovQ;
        end
        
        newq(zIx, bIx) = dot(zPi(zIx, :), expTerms) / (1.0 + r);
      end
    end % zIx
  end % bIx
  
  err = max(abs(q(:)  - newq(:) ));
  q = newq;
  if mod(iter,10)==0
    fprintf('Iter %d with error %f \n', iter, err);
  end
  iter = iter + 1;
end

figure;
subplot(1, 2, 1); plot(b, bPr(1:5:end, :)); hold on; plot(b, b, 'k--', 'LineWidth', 2); title('b''');
% subplot(1, 3, 2); plot(b, d(1:3:end, :)); title('d');
qb = q .* repmat(b, [zSz, 1]);
subplot(1, 2, 2); plot(b, qb(1:5:end, :)); title('q b');

rng(1989);
T = 50000;
[ simZ, simZix ] = simulateMarkov( T, z, zPi );
simB = bBar * ones([T, 1]);
simD = zeros([T, 1]);
simQ = qMin * ones([T, 1]);

for tIx = 1:T-1
  bHere = simB(tIx);
  zHere = simZ(tIx);
  zIx = simZix(tIx);
  laff = q(zIx, :) .* (b - (1.0 - delta) * bHere); %#ok<*PFBNS>
  laff(q(zIx, :) < qMin) = -1000;
  thold = kappa * bHere - prim(bHere) - zHere;
  if max(laff) < thold
    if bHere > 0.0
      simD(tIx) = 1;
      simB(tIx+1) = min(maxRecov, nu * simB(tIx));
    else
      simD(tIx) = 0;
      simB(tIx+1) = 0.0;
    end
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

          slope = (q(zIx, bPrIx+1) - q(zIx, bPrIx)) / bInterval;
          intercept = q(zIx, bPrIx) - slope * b(bPrIx);
          simQ(tIx+1) = intercept + slope * simB(tIx+1);
          break;
        end
      end
    end
  end
end

simH = prim(simB) + simZ';
simSp = kappa * (1.0 ./ simQ - 1.0);

figure;
subplot(2, 2, 1); histogram(simZ, zSz, 'Normalization', 'probability'); hold on;
histogram(simH, 50, 'Normalization', 'probability'); title('z and H');
subplot(2, 2, 2); histogram(simB, 50); title('B');
subplot(2, 2, 3); histogram(simD, 2, 'Normalization', 'probability'); title('Default');
subplot(2, 2, 4); histogram(simSp, 80, 'Normalization', 'probability'); xlim([0.0 0.1]); title('Spread');

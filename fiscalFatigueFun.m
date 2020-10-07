clear; clc;

bBar = 0.3;
rf = 0.04;

atB = [ 0.0, bBar, 1.5];
primSurp = [ -0.05, rf*bBar, 0.08 ];

fatigue = @(coefs, b) coefs(1) +  coefs(2) ./ (1.0 + exp( - coefs(3) * b));

ccf = fsolve( @(cc) primSurp - fatigue(cc, atB), [0.0 0.5 1.0]);

scatter(atB, primSurp);
hold on;
fplot( @(bb) fatigue(ccf, bb), [-0.75, 1.75]);


%% Euler Eq Stuff

numsteps = [100, 1000, 10000];
taus = [0.1, 0.01, 0.001];

y = zeros(3, max(numsteps));
y(:,1) = 1;
dydts = -y(:,1);


for i = 1:(numsteps(1) - 1)
    y(1, i + 1) = y(1, i) + taus(1) * dydts(1);
    dydts(1) = -y(1, i + 1);
end

for i = 1:(numsteps(2) - 1)
    y(2, i + 1) = y(2, i) + taus(2) * dydts(2);
    dydts(2) = -y(2, i + 1);
end

for i = 1:(numsteps(3) - 1)
    y(3, i + 1) = y(3, i) + taus(3) * dydts(3);
    dydts(3) = -y(3, i + 1);
end

times1 = linspace(0, taus(1)*numsteps(1), numsteps(1));
times2 = linspace(0, taus(2)*numsteps(2), numsteps(2));
times3 = linspace(0, taus(3)*numsteps(3), numsteps(3));

y_exact1 = exp(-times1);
y_exact2 = exp(-times2);
y_exact3 = exp(-times3);

figure(1)
hold on
plot(times1, y(1,1:numsteps(1)))
plot(times2, y(2,1:numsteps(2)))
plot(times3, y(3,1:numsteps(3)))
plot(times3, y_exact3)
legend('tau = 0.1', 'tau = 0.01', 'tau = 0.001', 'analytical')
hold off

MeanError = zeros(1,3);
MeanError(1) = mean(y(1,1:numsteps(1)) - y_exact1);
MeanError(2) = mean(y(2,1:numsteps(2)) - y_exact2);
MeanError(3) = mean(y(3,1:numsteps(3)) - y_exact3);

figure(2)
loglog(taus, MeanError, '-o')
xlabel('Time Step Size')
ylabel('Mean Error')
grid on

% trying midpoint method

numstep = 100;
tau = 0.1;

ymid = 1;
dydt = -exp(-tau/2);

for i = 1:(numstep - 1)
    ymid(i + 1) = ymid(i) + tau * dydt;
    dydt = -exp(-i * tau - tau/2);
end

times = linspace(0, tau*numstep, numstep);

figure(3)
hold on
plot(times, ymid)
plot(times1, y_exact1)
legend('midpoint method', 'analytical')
hold off

MeanError2 = mean(ymid - y_exact1);

Error = ymid - y_exact1;
percerror = (ymid - y_exact1) ./ y_exact1 .* 100;

figure(4)
plot(times, percerror)

figure(5)
plot(times, Error)

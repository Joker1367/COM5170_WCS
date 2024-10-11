clear all;
clc;

v_c = 3e8;
N = 1e8;
theta = -pi + 2 * pi * rand(1, N);
%% (a)
f_c = 2e9;
v = 20 * 1000 / 3600;
f_m = v / (v_c / f_c);

f_doppler = f_m * cos(theta);

figure(1)
grid on;
histogram(f_doppler, 200, 'Normalization', 'pdf');
xlabel('Doppler shift (Hz)');
ylabel('Probability');
title('pdf of f_D under v = 20km/hr and f_c = 2GHz');

figure(2)
grid on;
histogram(f_doppler, 200, 'Normalization', 'cdf');
xlabel('Doppler shift (Hz)');
ylabel('Probability');
title('cdf of f_D under v = 20km/hr and f_c = 2GHz');

%% (b)
f_c = 26e9;
v = 90 * 1000 / 3600;
f_m = v / (v_c / f_c);

f_doppler = f_m * cos(theta);

figure(3)
grid on;
histogram(f_doppler, 200, 'Normalization', 'pdf');
xlabel('Doppler shift (Hz)');
ylabel('Probability');
title('pdf of f_D under v = 90km/hr and f_c = 26GHz');

figure(4)
grid on;
histogram(f_doppler, 200, 'Normalization', 'cdf');
xlabel('Doppler shift (Hz)');
ylabel('Probability');
title('CDF of f_D under v = 90km/hr and f_c = 26GHz');

%% (c)
f_c = 2e9;
v = (20 + 70 * rand(1,N)) * 1000 / 3600;
f_m = v / (v_c / f_c);

f_doppler = f_m .* cos(theta);

figure(5)
grid on;
histogram(f_doppler, 200, 'Normalization', 'pdf');
xlabel('Doppler shift (Hz)');
ylabel('Probability');
title('pdf of f_D under v ~ U(20,90) km/hr and f_c = 2GHz');

figure(6)
grid on;
histogram(f_doppler, 200, 'Normalization', 'cdf');
xlabel('Doppler shift (Hz)');
ylabel('Probability');
title('pdf of f_D under v ~ U(20,90) km/hr and f_c = 2GHz');

%% (d)
f_c = 2e9;
v = 20 * 1000 / 3600;
f_m = v / (v_c / f_c);

f_doppler = f_m * cos(theta);

f_D = linspace(-f_m, f_m, 1e3);
pdf = (1 / pi) * (1 ./ sqrt(f_m * f_m - f_D .* f_D));
cdf = (1 / pi) * asin(f_D / f_m) + 0.5;

figure(7)
hold on;
grid on;
plot(f_D, pdf, 'b-', 'LineWidth', 2);
histogram(f_doppler, 200, 'Normalization', 'pdf');
xlabel('Doppler shift (Hz)');
ylabel('Probability');
legend('theoretical', 'simulated')
title('pdf of f_D under v = 20km/hr and f_c = 2GHz');

figure(8)
hold on;
grid on;
plot(f_D, cdf, 'b-', 'LineWidth', 2);
histogram(f_doppler, 200, 'Normalization', 'cdf');
xlabel('Doppler shift (Hz)');
ylabel('Probability');
legend('theoretical', 'simulated')
title('cdf of f_D under v = 20km/hr and f_c = 2GHz');
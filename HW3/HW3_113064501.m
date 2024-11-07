clear all;
clc;

%% Q1 Filtered Gaussian Noise Method

rng(0);
p = 1;

Omega_p = 1;
fmT = [0.01 0.1 0.5];
T = 1e6;

zeta = [];
sigma_2 = [];

duration = 300;
starting = 2001;
ending = starting + duration - 1;

for i = 1 : length(fmT)
    kp = 2 - cos(pi * fmT(i) / 2) - sqrt((2 - cos(pi * fmT(i) / 2))^2 - 1);
    zeta = [zeta kp];
    sigma_2 = [sigma_2 (Omega_p / 2) * (1 + kp) / (1 - kp)];
end

g_I = zeros(length(fmT), T);
g_Q = zeros(length(fmT), T);
samples = zeros(length(fmT), T);

for i = 1 : length(fmT)
    g_I(i, 1) = sqrt(Omega_p / 2) * randn;
    g_Q(i, 1) = sqrt(Omega_p / 2) * randn;
    samples(i, 1) = sqrt(g_I(i, 1) * g_I(i, 1) + g_Q(i, 1) * g_Q(i, 1));
    for t = 2 : T
        g_I(i, t) = zeta(i) * g_I(i, t - 1) + (1 - zeta(i)) * sqrt(sigma_2(i)) * randn;
        g_Q(i, t) = zeta(i) * g_Q(i, t - 1) + (1 - zeta(i)) * sqrt(sigma_2(i)) * randn;
        samples(i, t) = sqrt(g_I(i, t) * g_I(i, t) + g_Q(i, t) * g_Q(i, t));
    end
end

samples = 10*log10(samples);

figure(p)
p = p+1;
subplot(311)
plot([1:300], samples(1,starting:ending));
xlabel('Time,t/T');
ylabel('Envelope Level(dB)');
title(sprintf('Channel Output for fmT=0.01'));
subplot(312)
plot([1:300], samples(2,starting:ending));
xlabel('Time,t/T');
ylabel('Envelope Level(dB)');
title(sprintf('Channel Output for fmT=0.1'));
subplot(313)
plot([1:300], samples(3,starting:ending));
xlabel('Time,t/T');
ylabel('Envelope Level(dB)');
title(sprintf('Channel Output for fmT=0.5'));

figure(p)
p = p + 1;
hold on
for i = 1 : length(fmT)
    auto = xcorr(g_I(i,:) + 1i*g_Q(i,:), 10/fmT(i), 'normalized');
    auto = auto((length(auto) + 1) / 2:end);
    axis = 0:fmT(i):10;
    plot(axis, auto,'-');
end
xlabel('Time delay,f_m\tau');
ylabel('Autocorrelation \phiII(\tau)');
title('Channel Output autocorrelation ');
legend('f_mT=0.01','f_mT=0.1','f_mT=0.5');

%% Q2 Sum of Sinusoids method
fmT = [0.01, 0.1, 0.5];
M   = [8, 16];
T = 1e6;

normalized_delay = 10;
duration = 300;
starting = 2001;
ending = starting + duration - 1;

for i = 1 : length(fmT)
    g_I = zeros(length(M), T);
    g_Q = zeros(length(M), T);
    samples = zeros(length(M), T);
    for j = 1 : length(M)
        N = 4 * M(j) + 2;
        for t = 1 : T
            f_m_t = fmT(i) * t;
            g_I(j, t) = 2 * cos(2 * pi * f_m_t);
            for n = 1 : M(j)
                beta_n = pi * n / M(j);
                f_n_t = f_m_t * cos(2 * pi * n / N);
                g_I(j, t) = g_I(j, t) + 2*sqrt(2)*cos(beta_n)*cos(2*pi*f_n_t);
                g_Q(j, t) = g_Q(j, t) + 2*sqrt(2)*sin(beta_n)*cos(2*pi*f_n_t);
            end
            samples(j, t) = 10 * log10(sqrt(g_I(j, t) * g_I(j, t) + g_Q(j, t) * g_Q(j, t)));
        end
    end
    figure(p)
    p = p+1;
    subplot(211)
    plot([1:300], samples(1,starting:ending));
    xlabel('Time,t/T');
    ylabel('Envelope Level(dB)');
    title(sprintf('Channel Output for fmT=%.2f M=8', fmT(i)));
    subplot(212)
    plot([1:300], samples(2,starting:ending));
    xlabel('Time,t/T');
    ylabel('Envelope Level(dB)');
    title(sprintf('Channel Output for fmT=%.2f M=16', fmT(i)));

    figure(p)
    p = p + 1;
    hold on
    for j = 1 : length(M)
        auto = xcorr(g_I(j,:) + 1i*g_Q(j,:), 10/fmT(i), 'normalized');
        auto = auto((length(auto) + 1) / 2:end);
        axis = 0:fmT(i):10;
        plot(axis, auto,'-');
    end
    xlabel('Time delay,f_m\tau');
    ylabel('Autocorrelation \phiII(\tau)');
    title(sprintf('Channel Output autocorrelation for f_mT = %.2f', fmT(i)));
    legend('M = 8','M = 16');
end


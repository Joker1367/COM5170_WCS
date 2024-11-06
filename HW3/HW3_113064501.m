clear all;
clc;

%% 1-1

rng(0);
Omega_p = 1;
fmT = [0.01 0.1 0.5];
normalize_delay = 10;

zeta = [];
sigma_2 = [];

for i = 1 : length(fmT)
    kp = 2 - cos(pi * fmT(i) / 2) - sqrt((2 - cos(pi * fmT(i) / 2))^2 - 1);
    zeta = [zeta kp];
    sigma_2 = [sigma_2 (Omega_p / 2) * (1 + kp) / (1 - kp)];
end

g_I = zeros(length(fmT), normalize_delay / min(fmT) + 1);
g_Q = zeros(length(fmT), normalize_delay / min(fmT) + 1);
samples = zeros(length(fmT), normalize_delay / min(fmT) + 1);

for i = 1 : length(fmT)
    g_I(i, 1) = sqrt(Omega_p / 2) * randn;
    g_Q(i, 1) = sqrt(Omega_p / 2) * randn;
    samples(i, 1) = sqrt(g_I(i, 1) * g_I(i, 1) + g_Q(i, 1) * g_Q(i, 1));
    for t = 2 : normalize_delay / min(fmT) + 1
        g_I(i, t) = zeta(i) * g_I(i, t - 1) + (1 - zeta(i)) * sqrt(sigma_2(i)) * randn;
        g_Q(i, t) = zeta(i) * g_Q(i, t - 1) + (1 - zeta(i)) * sqrt(sigma_2(i)) * randn;
        samples(i, t) = sqrt(g_I(i, t) * g_I(i, t) + g_Q(i, t) * g_Q(i, t));
    end
end

samples = 10*log10(samples);

figure(1);
plot([1:300],samples(1,[1:300]));
xlabel('Time,t/T');
ylabel('Envelope Level(dB)');
title('Channel Output for fmT=0.01');

figure(2);
plot([1:300],samples(2,[1:300]));
xlabel('Time,t/T');
ylabel('Envelope Level(dB)');
title('Channel Output for fmT=0.1');

figure(3);
plot([1:300],samples(3,[1:300]));
xlabel('Time,t/T');
ylabel('Envelope Level(dB)');
title('Channel Output for fmT=0.5');

%% 1-2
figure(4)
hold on
for i = 1 : length(fmT)
    auto = xcorr(g_I(i,:), 10/fmT(i), 'normalized');
    auto_gI = auto((length(auto) + 1) / 2:end);
    axis = 0:fmT(i):10;
    plot(axis, auto_gI,'-');
end
xlabel('Time delay,fm\tau');
ylabel('Autocorrelation \phiII(\tau)');
title('Channel Output autocorrelation ');
legend('f_m_T=0.01','f_m_T=0.1','f_m_T=0.5');

%% 2-1
clear all
clc

fmT = [0.01, 0.1, 0.5];
M   = [8, 16];
normalize_delay = 10;

p = 5;

for i = 1 : length(M)
    g_I = zeros(length(fmT), normalize_delay / min(fmT) + 1);
    g_Q = zeros(length(fmT), normalize_delay / min(fmT) + 1);
    samples = zeros(length(fmT), normalize_delay / min(fmT) + 1);
    auto = zeros(length(fmT), normalize_delay / min(fmT) + 1);
    for j = 1 : length(fmT)
        N = 4 * M(i) + 2;
        for t = 1 : normalize_delay / min(fmT) + 1
            f_m_t = fmT(j) * t;
            g_I(j, t) = 2 * cos(2 * pi * f_m_t);
            for n = 1 : M(i)
                beta_n = pi * n / M(i);
                f_n_t = f_m_t * cos(2 * pi * n / N);
                g_I(j, t) = g_I(j, t) + 2*sqrt(2)*cos(beta_n)*cos(2*pi*f_n_t);
                g_Q(j, t) = g_Q(j, t) + 2*sqrt(2)*sin(beta_n)*cos(2*pi*f_n_t);
            end
            samples(j, t) = 10 * log10(sqrt(g_I(j, t) * g_I(j, t) + g_Q(j, t) * g_Q(j, t)));
        end
    end
    for j = 1 : length(fmT)
        tmp = xcorr(g_I(i,:) + 1i*g_Q(i,:), 10/fmT(i), 'normalized');
        auto(j,:) = tmp((length(auto) + 1) / 2:end);
    end
    figure(p)
    p = p+1;
    subplot(311)
    plot([1:300], samples(1,1:300));
    xlabel('Time,t/T');
    ylabel('Envelope Level(dB)');
    title(sprintf('Channel Output for fmT=0.01 M=%d', M(i)));
    subplot(312)
    plot([1:300], samples(2,1:300));
    xlabel('Time,t/T');
    ylabel('Envelope Level(dB)');
    title(sprintf('Channel Output for fmT=0.1 M=%d', M(i)));
    subplot(313)
    plot([1:300], samples(3,1:300));
    xlabel('Time,t/T');
    ylabel('Envelope Level(dB)');
    title(sprintf('Channel Output for fmT=0.1 M=%d', M(i)));


end


%% 2-2

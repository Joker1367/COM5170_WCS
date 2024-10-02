clear all;
clc;

%%
syms rho v;

block_prob = [0.01, 0.03, 0.05, 0.1];

result = zeros(41, 4);

%% Q1
for m = 1 : 20
    cur = 1;
    denominator = 1;
    for k = 1 : m
        cur = cur * (rho / k);
        denominator = denominator + cur;
    end

    log_factorial = 0;
    for k = 1 : m
        log_factorial = log_factorial + log(k);
    end

    for block_prob_idx = 1 : 4
        equation = log(block_prob(block_prob_idx)) == log(rho^m) - log_factorial - log(denominator);
        solution = vpasolve(equation);
        result(m, block_prob_idx) = double(max(real(solution)));
    end
end

for m = 200:220
    cur = 1;
    denominator = 1;
    for k = 1 : m
        cur = cur * (rho / k);
        denominator = denominator + cur;
    end
    
    log_factorial = 0;
    for k = 1 : m
        log_factorial = log_factorial + log(k);
    end

    for block_prob_idx = 1 : 4
        equation = log(block_prob(block_prob_idx)) == log(rho^m) - log_factorial - log(denominator);
        solution = vpasolve(equation);
        result(m - 199 + 20, block_prob_idx) = double(max(real(solution)));
    end
end

%% Q3
operator = [1, 2, 3];
traffic_load = zeros(3, 4);
trunking_efficiency = traffic_load;

for operator_idx = 1:length(operator)
    m = 600 / (5 * operator(operator_idx));
    cur = 1;
    denominator = 1;
    for k = 1 : m
        cur = cur * (rho / k);
        denominator = denominator + cur;
    end
    
    log_factorial = 0;
    for k = 1 : m
        log_factorial = log_factorial + log(k);
    end

    for block_prob_idx = 1 : length(block_prob)
        equation = log(block_prob(block_prob_idx)) == log(rho^m) - log_factorial - log(denominator);
        solution = vpasolve(equation);
        traffic_load(operator_idx, block_prob_idx) = double(max(real(solution)));
        trunking_efficiency(operator_idx, block_prob_idx) = traffic_load(operator_idx, block_prob_idx) / m;
    end
end
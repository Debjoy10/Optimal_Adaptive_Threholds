% CPS Group-2: Fixed threshold attack detector for dynamical systems
%Gives Optimal Delay for fixed threshold algorithm
function [optimal_delay] = fixed_threshold(~)
    timeWindow = 15;
    delay = 0;
    optimal_loss = Inf;
    C = 0.1;
    D_table = readtable('flies/delay_th=0.05.csv');
    D_array = table2array(D_table);
    P = zeros(timeWindow,timeWindow);
    while delay < timeWindow
        k_a = 1;
        P_dash = 0;
        L_dash = inf;
        while k_a < (timeWindow - delay)
            P(delay+1,k_a) = D_array(k_a+delay,k_a);
            if P(delay+1,k_a) > P_dash
                P_dash = P(delay+1,k_a);
                L_dash = P_dash + C*false_positive_rate(delay)*timeWindow;
            end
            k_a = k_a + 1;
        end
        if L_dash < optimal_loss
            optimal_loss = L_dash;
            optimal_delay = delay;
        end
        delay = delay + 1;
    end
end

% Simple function for FP rate
function FP = false_positive_rate(del)
    FP = 1 - del/15;
end

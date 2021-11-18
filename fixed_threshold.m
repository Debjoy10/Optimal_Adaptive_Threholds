% CPS Group-2: Fixed threshold attack detector for dynamical systems
%Gives Optimal Delay for fixed threshold algorithm
function [optimal_threshold] = fixed_threshold(~)
    timeWindow = 15;
    delay = 0;
    optimal_loss = Inf;
    C = 0.1;
    D_table = readtable('files/delay_th=0.05.csv');
    D_array = table2array(D_table);
    load('files/delay_to_thresh', 'd_to_th')
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
    disp(optimal_delay)
    threshold = d_to_th(optimal_delay);
    disp(threshold)
    threshold_array = ones(1,timeWindow);
    optimal_threshold = threshold.*threshold_array;
    writematrix(optimal_threshold,'files/fixed_threshol.csv')
end

% Simple function for FP rate
function FP = false_positive_rate(del)
    FP = 1 - del/15;
end

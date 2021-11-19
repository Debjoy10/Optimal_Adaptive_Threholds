function [optimal_threshold] = fixed_threshold(~)
    timeWindow = 15;
    optimal_loss = Inf;
    C = 0.1;
    load('files/thresh_to_damage', 'th_to_damage')
    load('files/delay_to_thresh', 'd_to_th')
    load('files/all_delays', 'all_delays');
    P = zeros(timeWindow,timeWindow);
    for delay = all_delays
        k_a = 1;
        P_dash = 0;
        L_dash = inf;
        damage = th_to_damage(d_to_th(delay));
        while k_a < (timeWindow - delay)
            P(delay+1,k_a) = damage(k_a+delay,k_a);
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
    end
    disp(optimal_delay)
    threshold = d_to_th(optimal_delay);
    disp(threshold)
    threshold_array = ones(1,timeWindow);
    optimal_threshold = threshold.*threshold_array;
    writematrix(optimal_threshold,'files/fixed_threshold.csv')
end

% Simple function for FP rate
function FP = false_positive_rate(del)
    FP = 1 - del/15;
end
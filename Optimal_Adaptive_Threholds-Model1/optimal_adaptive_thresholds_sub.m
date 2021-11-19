% CPS Group-2: Adaptive attack detector for dynamical systems
% Runs the optimal adaptive threholds algorithm

function [TCP_opt, optimal_delays, optimal_thresholds] = optimal_adaptive_thresholds_sub(~)
    timeWindow = 15;
    % Set of all damage values and delay to threshold mappings
    load('files/all_damages', 'all_damages')
    load('files/delay_to_thresh', 'd_to_th')
    
    % Run the MinimumCostThresholds algorithm over all possible damage values
    TCP_opt = Inf;
    for P = all_damages
        [TC_P, optimal_delays_P] = MinimumCostThresholds(P, timeWindow);
        if TC_P < TCP_opt
            TCP_opt = TC_P;
            optimal_delays = optimal_delays_P;
        end
    end
    
    % Get the optimal thresholds from delay - threshold mapping
    optimal_thresholds = zeros(timeWindow+1, 1);
    for idx = 1:timeWindow+1
        optimal_thresholds(idx) = d_to_th(optimal_delays(idx));
    end
end

function [cost, delay_opt] = MinimumCostThresholds(P, T)
    timeWindow = T;
    C = 10;  % FP cost
    Cd = 0; % Threshold Change Cost
    
    % Load all mapping files
    load('files/thresh_to_damage', 'th_to_damage')
    load('files/delay_to_thresh', 'd_to_th')
    load('files/all_delays', 'all_delays');
    
    % To store cost and optimal delays
    COST = double(1000*ones(timeWindow+1, timeWindow+1, timeWindow+1));
    del_opt = double(zeros(timeWindow+1, timeWindow+1, timeWindow+1));
    for m = 0:timeWindow-1
        for del = all_delays
            COST(timeWindow+1, m+1, del+1) = 0;
        end
    end
    
    for n = timeWindow:-1:1
        for m = 0:n-1
            for del_n_1 = all_delays
                S = double(1000*ones(max(all_delays), 1));
                for del_n = all_delays 
                    damage = th_to_damage(d_to_th(del_n));
                    if del_n > m
                        S(del_n+1) = COST(n+1,m+2,del_n+1) + C*false_positive_rate(del_n);
                    elseif damage(n, n-m) <= P
                        S(del_n+1) = COST(n+1,del_n+1,del_n+1) + C*false_positive_rate(del_n);
                    else
                        S(del_n+1) = 1000;
                    end
%                     disp(del_n)
%                     disp(d_to_th(del_n))
%                     disp(false_positive_rate(del_n))
%                     disp(S(del_n+1))
%                     disp("---")
                    % Add threshold change cost (if applicable)
                    if (del_n_1 ~= del_n) && n > 1
                        S(del_n+1) = S(del_n+1) + Cd;
                    end
                end
                
                % Find Optimal Delays and Cost
                [argvalue, argmin] = min(S);
                del_opt(n, m+1, del_n_1+1) = argmin-1;
                COST(n, m+1, del_n_1+1) = argvalue;
            end
        end
    end
    writematrix(COST, 'files/COST.csv');
    m = 0;
    delay_opt = zeros(timeWindow+1, 1);
    delay_opt(1) = randperm(timeWindow,1); % Arbitrary Value
    for n = 1:timeWindow
        delay_opt(1+n) = del_opt(n, m+1, delay_opt(n)+1);
        m = min(m+1, delay_opt(1+n));
    end
    cost = COST(1, 1, delay_opt(1));
end

% Simple function for FP rate
function FP = false_positive_rate(del)
    FP = double(1 - (del/15));
end
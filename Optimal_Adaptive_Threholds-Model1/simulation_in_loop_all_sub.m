% CPS Group-2: Adaptive attack detector for dynamical systems
% Runs simulation and stores delay/damage values
% Simulation for different threhold values
% Constant Kalmann and LQR gains
% Both Detectors Supported

% System: motor position control
J = 3.2284E-6;
b = 3.5077E-6;
KK = 0.0274;
RR = 4;
LL = 2.75E-6;

Ac = [0 1 0; 0 -b/J KK/J; 0 -KK/LL -RR/LL];
Bc = [0; 0; 1/LL];
Cc = [1 0 0];
Dc = [0];

x0 = [5; 0; -1000]; % initial state
Ts= 0.1;
ref = 0; % ref??
sys_c = ss(Ac,Bc,Cc,Dc);
sys = c2d(sys_c,Ts);

A = sys.a;
B = sys.b;
C = sys.c;
D = sys.d;


Q = [ 1 0 0 ;
      0 0.1 0;
      0 0 1;
    ];

R = 0.1;
R = R*Ts*Ts;
Q = Q*Ts*Ts;
N = zeros(size(Ac,2),size(Bc,2));
[K] = dlqr(A,B,Q,R,N);
% K = [16.0302    5.6622]; % LQR gain

QN = 1;
RN = 10;
sys_gain = ss(A- B*K,B,C,D,Ts);
[kalmf,L,P] = kalman(sys,QN,RN);

% A= [1.0000    0.1000;0    1.0000];
% B= [0.0050;0.1000];
% C= [1 0];
% D= [0];
size_x = [3 1];
size_y = [1 1]; 

timeWindow = 25;
sensorAttack = 0.2;
actuatorAttack = 0.2;
safex = [4;20;20];
depth = 0.1;
cusum_true = true;
cusum_cost_mat = [1]; %In case Y is also a vector, then we would require to normalize it

th_all = [1.8:0.2:3 3.5:0.5:6 8:1:10];

% sensorRange = [-4;4];
% actuatorRange = [-5;5]; 
% safex = [-4,-20,0;4,20,20];

% sensorRange = 4;
% actuatorRange = 5;

for th = th_all
    safex = [4;20;20]; % 3*1
    depth = 0.1;
    timeWindow = 25;
    sensorAttack = 0.5;%2;
    actuatorAttack = 0.5;%2;

    delay = 100*ones(timeWindow,timeWindow);
    damage = 100*ones(timeWindow,timeWindow);
    for k_a=2:timeWindow %attack start
        for k_e=1:k_a-1 %attack end
            delay(k_e,k_a) = 100; %a large value indicating attack not possible in this sub-window
        end
    end

    % initialize
    x_a = x0;
    xhat_a = x0;
    u_a = -K*xhat_a;

    x = x0;
    xhat = x0;
    u = -K*xhat;

    for k_a=1:timeWindow % attack start
        for k_e=k_a:timeWindow % attack end
            d = 0;
            p = 0;
            S_p = zeros(size_x);
            S_n = zeros(size_x);
            detected = false;
            for i=1:timeWindow
                % Non-attack scenario
                x = A*x + B*u; % state updattion in plant side
                y = C*x; %sensor output % sensor measurement in plant side
                r = y - C*xhat; % residue computation in controller side
                xhat = A*xhat + B*u + L*r; % state estimation in controller side
                u = - K*xhat; % control signal computation in controller side
                u_list(:,i) = u;
                % attack scenario
                if i>=k_a && i<=k_e % attack window
                    x_a = A*x_a + B*(u_a + actuatorAttack);
                else
                    x_a = A*x_a + B*u_a;
                end
                
                y_a = C*x_a; %sensor output
                if i>=k_a && i<=k_e % attack window
                    y_a = y_a + sensorAttack;
                end
                r_a = y_a - C*xhat_a; % residue
                xhat_a = A*xhat_a + B*u_a + L*r_a;
                u_a = - K*xhat_a;

                if cusum_true 
                      %S_p = max([0, S_p + r_a]); %If residue is +ve, it is added to S_p, if it was -ve, then S_p drops to a min 0
                      %S_n = min([0, S_n + r_a]); %If residue is -ve, it is added to S_n, making it more -ve. If it is +ve, S_n increases, moving towards 0 
                      for j = 1:size_y(1)
                          S_p(j) = max(0,S_p(j) + r_a(j));
                          S_n(j) = min(0,S_n(j) + r_a(j));
                      end
                      S_p_single = cusum_cost_mat*abs(S_p);
                      S_n_single = cusum_cost_mat*abs(S_n);
                      if(max(S_p_single,S_n_single)<th) 
                          if i >= k_a
                              %We also need to count the damage in this case as
                              %it will be plotted.
                              d = d + 1;
                              p = p + norm(x-x_a,inf);
                          end
                      else
                        %Attack is detected, need to break the loop.
                        detected = true;
                        break;
                      end
                else
                    if norm(r_a,inf)<th
                       if i >= k_a% you can use any norm here. But CUSUM is suggested as it is given in the paper
                           d = d + 1;
                           p = p + norm(x-x_a,inf);
                       end
                    else
                        detected = true;
                        break;
                    end
                end
            end
            if detected == true
                % Store Delay and Payoff
                delay(k_e,k_a) = d;
                damage(k_e,k_a) = p;
            end
        end
    end
    
    % Save Threshold Values in CSVs
    fout = sprintf('files/delay_th=%.2f.csv', th);
    writematrix(delay, fout);
    fout = sprintf('files/damage_th=%.2f.csv', th);
    writematrix(damage, fout);
end

% Extract Delay to Threshold Mapping
d_to_th = containers.Map('KeyType','double','ValueType','double');

% Store delay - th = max threshold that can lead to delay
for th = sort(th_all)
    fcsv = sprintf('files/delay_th=%.2f.csv', th);
    delay = readmatrix(fcsv);
    unique_delays = unique(delay);
    for d = unique_delays.'
        d_int = typecast(d, 'uint32');
        if d ~= 100
            d_to_th(d) = th;
        end
    end
end
all_delays = cell2mat(keys(d_to_th));

% Store Damage Values
P_all = -100;
th_to_damage  = containers.Map('KeyType', 'double', 'ValueType', 'any');
for th = sort(th_all)
    fcsv = sprintf('files/damage_th=%.2f.csv', th);
    P = readmatrix(fcsv);
    th_to_damage(th) = P;
    if P_all == -100
        P_all = P;
    else
        P_all = cat(1, P_all, P);
    end 
end
all_damages = unique(P_all).';

% Saving as Files
% 1. Delay to max-threshold that achieves that delay
% 2. Threshold to Damage Payoff Matrix
% 3. All Delay Values
% 4. All Damage Values
save('files/delay_to_thresh', 'd_to_th'); % load('files/delay_to_thresh', 'd_to_th')
save('files/thresh_to_damage', 'th_to_damage'); % load('files/thresh_to_damage', 'th_to_damage')
save('files/all_delays', 'all_delays');
save('files/all_damages', 'all_damages');


keys(d_to_th)
values(d_to_th)
Aa = [0.4450 -0.0458;
      1.2939 0.4402];
Bb = [0.0550;
      4.5607];
Cc = [0 1];
Dd = 0;
esp = ss(Aa,Bb,Cc,Dd);
sys_ss = esp;
x0 = [10;1]; % initial state
Ts = 0.04                                   ;
ref = 90; % reference to be tracked
p = 50000;
Q = 1*p*(C'*C);
R = 0.01;
[K,S,CLP] = dlqr(A,B,Q,R)
[xx,kk,ll]=icare(A,B,Q,R)
sys_cl = ss(A-B*K,B,C,D,Ts)
% [A,B,C,D]=ssdata(sys_cl);
x = x0;
u = -K*x;
for i=1:time
    x = A*x + B*u;
    u = -K*x;
    plot_vectorx1(i) = x(1);
    plot_vectorx2(i) = x(2);
    plot_vectoru(i) = u;
end
disp('Closed loop plant response in LQR');
clf
figure(1)
hold on
plot(plot_vectorx1)
plot(plot_vectorx2)
plot(plot_vectoru)
legend({'x1','x2','u'});
grid on
hold off
step(sys_cl);
initial(sys_cl,x0)
t=0:Ts:time;
reference = 1*ones(size(t));
lsim(sys_cl,reference,t,x0);
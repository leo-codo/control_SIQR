%Simulation SIQR model

%parameters

t = 180;
alpha = 0.2;
gamma  = 0.1;
mu = 0.02;
rho = 0.3;
epsilon = 0.1;
eta = 0.2;
Delta = 0.2;
v = 0.05;

% R_0

R0 = (Delta*alpha)/((mu + v)*(gamma + mu + eta))

%Equations

S = [9];
I = [1];
Q = [0];
R = [0];

N = [S(1) + I(1) + Q(1) + R(1)];

for i = 1:t,
  S(i+1) = S(i) + Delta - alpha*S(i)*I(i) - mu*S(i) - v*S(i);
  I(i+1) = I(i) + alpha*S(i)*I(i) - (gamma + mu + eta)*I(i);
  Q(i+1) = Q(i) + (eta - epsilon)*I(i) - (rho + mu)*Q(i);
  R(i+1) = R(i) + gamma*I(i) + rho*Q(i) - mu*R(i);
end;

%print results

S(t)
I(t)
Q(t)
R(t)

%Plots

figure (1);
hold on
for j=1:t
    plot(j,S(j),'*','color','r');
end
xlim([0,30])
title ("Susceptible population over time");
xlabel ("Time in days");
ylabel ("S(t)");
hold off

figure (2);
hold on
for j=1:t
    plot(j,I(j),'*','color','r');
end
xlim([0,20])
title ("Infectious population over time");
xlabel ("Time in days");
ylabel ("I(t)");
hold off

figure (3);
hold on
for j=1:t
    plot(j,Q(j),'*','color','r');
end
xlim([0,30])
title ("Quarantine population over time");
xlabel ("Time in days");
ylabel ("Q(t)");
hold off

figure (4);
hold on
for j=1:t
    plot(j,R(j),'*','color','r');
end
xlim([0,100])
title ("Removed population over time");
xlabel ("Time in days");
ylabel ("R(t)");
hold off

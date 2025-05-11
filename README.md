# control_SIQR

% Parameters
N = 100;
delta1 = 0.2 %0.00645/365;
alfa = 0.2 %0.2599;
pho =0.3 % 0.0619;
eta = 0.2% 0.001;
gamma = 0.1 %1/14.25;
mu = 0.02 %0.02056/365;
Death = 0.1%0.034;
v = 0.05;
R0 = (delta1*alfa)/((mu + v)*(gamma + mu + eta))
r1 = delta1/(mu+v)
% Initial conditions
S0 = 9;
I0 = 1;
Q0 = 0;
R0 = 0;

% Time vector
t = linspace(0, 365, 100);

% ODE solver
y0 = [S0, I0, Q0, R0];
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);
[t, solution] = ode45(@(t, y) SIQR_model(t, y, N, delta1, alfa, eta, pho, gamma, mu, Death, v), t, y0, options);

% Extracting solution
S = solution(:, 1);
I = solution(:, 2);
Q = solution(:, 3);
R = solution(:, 4);

% Plotting
figure;

subplot(2, 2, 1);
plot(t, S, 'r', 'LineWidth', 1);
title('Susceptible population over time');
xlabel('Time (days)');
ylabel('S(t)');
grid on;
subplot(2, 2, 2);
plot(t, I, 'b', 'LineWidth', 1);
title('Infectious population over time');
xlabel('Time (days)');
ylabel('I(t)');
grid on;
subplot(2, 2, 3);
plot(t, Q, 'g', 'LineWidth', 1);
title('Quarantine population over time');
xlabel('Time (days)');
ylabel('Q(t)');
grid on;
subplot(2, 2, 4);
plot(t, R, 'm', 'LineWidth', 1);
title('Recovered population over time');
xlabel('Time (days)');
ylabel('R(t)');
grid on;
% Print results
fprintf('S(t) at t=20 days: %f\n', S(30));
fprintf('I(t) at t=20 days: %f\n', I(30));
fprintf('Q(t) at t=20 days: %f\n', Q(end));
fprintf('R(t) at t=20 days: %f\n', R(end));

function dydt = SIQR_model(t, y,N, delta1, alfa, eta, pho, gamma, mu, Death, v)
    S = y(1);
    I = y(2);
    Q = y(3);
    R = y(4);
    
    dS_dt = delta1 - alfa * S * I - mu * S - v * S;
    dI_dt = alfa * S * I - (gamma + mu + eta) * I;
    dQ_dt = (eta - Death) * I - (pho + mu) * Q;
    dR_dt = pho * Q + gamma * I - mu * R;
    
    dydt = [dS_dt; dI_dt; dQ_dt; dR_dt];
end

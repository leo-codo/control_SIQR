% Definindo os valores das matrizes e parâmetros
mu = 0.02% 5.6329e-05;
v = 0.15;
alpha = 0.2%0.2599; %para R0=3.1;
delta = 0.2%1.7671e-05;
gamma = 0.1%0.0702;
eta = 0.2% 0.001;
epsilon = 0.1% 0.034;
rho = 0.3%0.0619;
lambda = 0.5;

% Definindo as matrizes A, B, R e Q
A = [-mu, 0, 0, 0;
      0, -(gamma + mu + eta), 0, 0;
      0, (eta - epsilon), -(rho + mu), 0;
      0, gamma, rho, -mu];

B = [1, 0;
      0, 1;
      0, 0;
      0, 0];

R = [2, 0;
      0, 2];

Q = [1, 0, 0, 0;
      0, 1, 0, 0;
      0, 0, 0, 0;
      0, 0, 0, 0];
%tempo final
T = 60;
% Resolvendo numericamente a equação de Riccati
tspan = [0 T]; % intervalo de tempo de 0 a 10
P0 = zeros(16, 1); % condição inicial para P
[t, P] = ode45(@(t, P) equacaoricatti(t, P, A, B, R, Q), tspan, P0);

% Calculando a norma de cada matriz P ao longo do tempo
norm_P = zeros(size(P, 1), 1);
for i = 1:size(P, 1)
    norm_P(i) = norm(reshape(P(i, :), [4, 4]));
end

% Plotando a norma da matriz P ao longo do tempo
figure;
plot(t, norm_P);
xlabel('Tempo');
ylabel('Norma de P');
title('Evolution of state variables over time');
grid on;



% Resolvendo numericamente a equação diferencial original
[t, x] = ode45(@(t, x) system_equation(t, x, A, B, v, alpha, delta), tspan, ones(4, 1));

% Extrair as variáveis de estado S, I, Q, R do vetor x
S = x(:, 1);
I = x(:, 2);
Q = x(:, 3);
R = x(:, 4);

% Plotando as variáveis de estado ao longo do tempo
figure;
plot(t, x);
xlabel('Tempo');
ylabel('State Variables');
legend('S', 'I', 'Q', 'R');
title('Evolution of state variables over time');
grid on;

% Print results
fprintf('S(t) at t=60 days: %f\n', S(end));
fprintf('I(t) at t=60 days: %f\n', I(end));
fprintf('Q(t) at t=60 days: %f\n', Q(end));
fprintf('R(t) at t=60 days: %f\n', R(end));




% Função para equação diferencial de Riccati
function dPdt = equacaoricatti(t, P, A, B, R, Q)
    P = reshape(P, [4, 4]);
    dPdt = -(A.' * P + P * A - P * B * inv(R) * B.' * P + Q);
    dPdt = dPdt(:);
end

% Função para equação diferencial original
function dxdt = system_equation(t, x, A, B, v, alpha, delta)
    mu = 0.02; % Definindo mu novamente aqui
    dxdt = (A * x) + B * [-v * x(1) - (alpha * delta / (mu + v)) * x(2);
                            (alpha * delta / (mu + v)) * x(2)];
end

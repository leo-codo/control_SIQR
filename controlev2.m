function main()
    % Parâmetros do sistema
    mu = 0.02;%0.00563;
%v = 10;
alpha = 0.2; %0.1606; %para R0=3.1;
Delta = 0.2;% 1.7671e-05;
gamma =0.1;% 0.0702;
eta = 0.2;
epsilon = 0.1;%0.034;
rho = 0.3;%0.0773;

    % Condições iniciais
    S0 = 9;
    I0 = 1;
    Q0 = 0;
    R0 = 0;

    % Tempo final
    T = 30;

    % Parâmetros do custo
    a1 = 1/2;
    a2 = 1/2;

    % Resolvendo as EDOs
    tspan = [0 T];
    y0 = [S0; I0; Q0; R0];
    [t, y] = ode45(@(t,y) ode_system(t, y, Delta, alpha, gamma, eta, epsilon, rho, mu, @control_function), tspan, y0);

    % Calculando o controle ótimo
    options = optimset('Display', 'off');
    v_opt = fminbnd(@(v) cost_function(t, y, v, a1, a2), -1, 1, options);

    % Exibindo o resultado
    disp(['Controle ótimo encontrado: v = ' num2str(v_opt)]);

    % Plotando as soluções
    figure;
    plot(t, y(:,1), 'b', t, y(:,2), 'r', t, y(:,3), 'g', t, y(:,4), 'm');
    legend('S', 'I', 'Q', 'R');
    xlabel('Time');
    ylabel('Population');
    title('System Dynamics');
    grid on
end

function dydt = ode_system(t, y, Delta, alpha, gamma, eta, epsilon, rho, mu, control_function)
    S = y(1);
    I = y(2);
    Q = y(3);
    R = y(4);

    v = control_function(t, y); % Obtendo o valor do controle
    dSdt = -(Delta*alpha/(mu+v) )* I - mu * S - v * S;
    dIdt = (Delta*alpha/(mu+v) )* I - (gamma + mu + eta) * I;
    dQdt = (eta - epsilon) * I - (rho + mu) * Q;
    dRdt = gamma * I + rho * Q - mu * R;

    dydt = [dSdt; dIdt; dQdt; dRdt];
end

function v = control_function(t, y)
    % Aqui você pode implementar a lógica para o controle ótimo
    % Por enquanto, retornamos um valor constante como exemplo
    v = 0.05;
                          
end

function cost = cost_function(t, y, v, a1, a2)
    S = y(:,1);
    I = y(:,2);

    cost_integrand = a1 * S.^2 + a2 * I.^2 + 0.5 * v.^2;
    cost = trapz(t, cost_integrand);
end
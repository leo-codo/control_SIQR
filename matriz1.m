% Definição da matriz A
mu = 0.00563;
v = 0.05;
alpha = 0.08;
Delta = 1;
gamma = 1/14;
eta = 0.2;
epsilon = 0.034;
rho = 0.3;

A = [-mu , 0, 0, 0;
    0, -(gamma + mu + eta), 0, 0;
    0, eta - epsilon, -(rho + mu), 0;
    0, gamma, rho, -mu];
% Definição da matriz B corrigida
B = [1 0;0 1;0 0;0 0];

% Verificação da controlabilidade
Wc = [B, A*B, A^2*B, A^3*B];
rank_Wc = rank(Wc);

% Verificação do posto
if rank_Wc == 4
    disp('O sistema é completamente controlável.');
else
    disp('O sistema pode não ser completamente controlável.');
end

% Mostrar o posto da matriz Wc
disp(['Posto da matriz Wc: ' num2str(rank_Wc)]);

    % Geração dos valores singulares
singular_values = svd(Wc);

% Índices dos valores singulares maiores que zero
nonzero_indices = find(singular_values > 0);

% Plotagem dos valores singulares maiores que zero
figure;
stem(nonzero_indices, singular_values(nonzero_indices), 'LineWidth', 2);
xlabel('Columns');
ylabel('Singular values');
title('Singular values ​​of the Wc controllability matrix');
grid on;

    % Mostrar o posto da matriz Wc
    disp(['Posto da matriz Wc: ' num2str(rank_Wc)]);

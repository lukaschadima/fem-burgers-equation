%% Příklad 1: Burgersova rovnice s přesným řešením (ověření) na [-1,1]
% Rovnice: u_t + u*u_x = nu * u_xx
% Přesné řešení dle (1.15)

% --- Parametry ---
nu = 0.01;              % Viskozita
c = 5;                  % Konstanta c > 1 pro hladkost
dom = [-1, 1];          % Doména
T = 1;                  % Konečný čas
tic;
% --- Definice PDE operátoru ---
pdefun = @(t, x, u) nu * diff(u, 2) - diff(u.^2 / 2);

% --- Okrajové podmínky (nulové Dirichletovy) ---
bc.left = @(u) u;       % u(-1) = 0
bc.right = @(u) u;      % u(1) = 0

% --- Počáteční podmínka ---
x = chebfun('x', dom);
u0 = (2 * nu * pi * sin(pi * x)) ./ (c + cos(pi * x));

% --- Časový vektor ---
t = linspace(0, T, 30);

% --- Spuštění pde15s s vysokou přesností ---
opts = pdeset('Eps', 1e-8);
[t_sol, u_sol] = pde15s(pdefun, t, u0, bc, opts);

% --- Výpočet přesného řešení ---
exact_solution = @(x,t) (2 * nu * pi * exp(-nu * pi^2 * t) .* sin(pi * x)) ./ ...
                        (c + exp(-nu * pi^2 * t) .* cos(pi * x));
toc;
% --- Vizualizace a chyba ---
figure;
subplot(1,2,1)
plot(u_sol(:,end), 'LineWidth', 2), hold on
plot(exact_solution(x, T), '--r', 'LineWidth', 2)
title('Porovnání řešení v čase T')
legend('Chebfun (pde15s)', 'Přesné', 'Location', 'best')
xlabel('x'), ylabel('u(x,T)'), grid on

subplot(1,2,2)
err = u_sol(:,end) - exact_solution(x, T);
plot(err, 'k', 'LineWidth', 2)
title(sprintf('Chyba (max = %.2e)', norm(err, inf)))
xlabel('x'), ylabel('Chyba'), grid on
%% Příklad 2: Burgersova rovnice se zdrojovým členem na [-1,1]
% Rovnice: u_t + u*u_x = nu * u_xx + f(x,t)

nu = 0.05;   omega = 2*pi;   dom = [-1, 1];   T = 10.67;
tic;
% --- Přesné řešení a zdroj ---
u_exact = @(x,t) sin(pi * x) .* sin(omega * t);
f_source = @(x, t) omega * sin(pi * x) .* cos(omega * t) + ...
                   (pi/2) * sin(2 * pi * x) .* sin(omega * t).^2 + ...
                   nu * pi^2 * sin(pi * x) .* sin(omega * t);

% --- PDE operátor včetně zdroje ---
pdefun = @(t, x, u) nu * diff(u, 2) - diff(u.^2 / 2) + f_source(x, t);

% --- Okrajové a počáteční podmínky ---
bc.left = @(u) u;   bc.right = @(u) u;
x = chebfun('x', dom);
u0 = 0 * x;          % nulová počáteční podmínka

t = linspace(0, T, 30);
opts = pdeset('Eps', 1e-14);
[t_sol, u_sol] = pde15s(pdefun, t, u0, bc, opts);
toc;
% --- Srovnání ---
u_cheb_T = u_sol(:,end);
u_ex_T = u_exact(x, T);

figure;
subplot(1,2,1)
plot(u_cheb_T, 'LineWidth', 2), hold on
plot(u_ex_T, '--r', 'LineWidth', 2)
title(sprintf('Porovnání v t = %.2f', T))
xlabel('x'), ylabel('u(x,T)')
legend('Chebfun (pde15s)', 'Přesné'), grid on

subplot(1,2,2)
err = u_cheb_T - u_ex_T;
plot(err, 'k', 'LineWidth', 2)
title(sprintf('Chyba (max = %.2e)', norm(err, inf)))
xlabel('x'), ylabel('Chyba'), grid on
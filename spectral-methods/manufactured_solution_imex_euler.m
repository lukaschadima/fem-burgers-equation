%% IMEX Euler pro Burgersovu rovnici se zdrojovým členem – Chebfun diskretizace na [-1,1]
% Rovnice: u_t + u*u_x = nu * u_xx + f(x,t)
% Přesné řešení: u(x,t) = sin(pi*x) * sin(omega*t)
clear; close all;

nu = 0.05;
omega = 2*pi;
dom = [-1, 1];
T = 1.67;
N = 64;                     % počet Čebyševových bodů
dt = 1e-2;                  % časový krok

% --- 1. Čebyševovy body a derivační matice ---
cheb_pts = chebpts(N, dom);
D1_raw = diffmat(N, 1, dom);
D2_raw = diffmat(N, 2, dom);

% --- 2. Vyřazení krajových bodů (okrajové podmínky u(-1)=u(1)=0) ---
idx = 2:N-1;
D1 = D1_raw(idx, idx);
D2 = D2_raw(idx, idx);
x_inner = cheb_pts(idx);    % vnitřní body (sloupcový vektor)

% --- 3. Zdrojový člen (funkce prostoru a času) ---
f_source = @(x, t) omega * sin(pi * x) .* cos(omega * t) + ...
    (pi/2) * sin(2 * pi * x) .* sin(omega * t).^2 + ...
    nu * pi^2 * sin(pi * x) .* sin(omega * t);

% --- 4. Počáteční podmínka (nulová) ---
u0_func = @(x) 0 * x;
u = u0_func(x_inner);       % vnitřní body

% --- 5. IMEX matice (pouze difúze implicitně) ---
I = eye(length(idx));
A = I - nu * dt * D2;       % (I - nu*dt*D2)

N_steps = ceil(T / dt);
t = 0;

for step = 1:N_steps
    % explicitní členy: advekce + zdroj v čase t
    u_x = D1 * u;
    conv = u .* u_x;                 % u * u_x
    f_val = f_source(x_inner, t);    % zdrojový člen v aktuálním čase

    rhs = u - dt * conv + dt * f_val;
    u = A \ rhs;
    t = t + dt;
end

% --- 6. Chyba vůči přesnému řešení ---
u_exact = @(x,t) sin(pi * x) .* sin(omega * t);
u_ref = u_exact(x_inner, T);
err = u - u_ref;
max_err = norm(err, inf);
fprintf('Max. chyba v t = %.2f:  %.3e\n', T, max_err);

% --- 7. Vykreslení ---
u_full = zeros(N, 1);
u_full(idx) = u;
u_cheb = chebfun(u_full, dom);
u_exact_cheb = chebfun(@(x) u_exact(x, T), dom);
err_cheb = u_cheb - u_exact_cheb;

figure;
subplot(1,2,1)
plot(u_cheb, 'b-', 'LineWidth', 1.5); hold on;
plot(u_exact_cheb, 'r--', 'LineWidth', 1.5);
xlabel('x'), ylabel('u(x,T)')
legend('IMEX Euler','Přesné','Location','best')
title(sprintf('t = %.2f (se zdrojem)', T)), grid on

subplot(1,2,2)
plot(err_cheb, 'k-', 'LineWidth', 1.5);
xlabel('x'), ylabel('Chyba')
title(sprintf('Chyba (max = %.2e)', max_err)), grid on
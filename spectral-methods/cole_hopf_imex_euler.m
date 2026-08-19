%% IMEX Euler pro Burgersovu rovnici – Chebfun diskretizace na [-1,1]
clear; close all;

nu = 0.01; c = 2; dom = [-1,1]; T = 5; N = 64; dt = 1e-1;

% ----- 1. Čebyševovy body a derivační matice -----
cheb_pts = chebpts(N, dom);             % body na [-1,1]
D1_raw = diffmat(N, 1, dom);            % 1. derivace
D2_raw = diffmat(N, 2, dom);            % 2. derivace

% ----- 2. Vyřazení krajových bodů (Dirichlet OK) -----
idx = 2:N-1;
D1 = D1_raw(idx, idx);
D2 = D2_raw(idx, idx);

% ----- 3. Počáteční podmínka na vnitřních bodech -----
u0_func = @(x) (2*nu*pi*sin(pi*x)) ./ (c + cos(pi*x));
u = u0_func(cheb_pts(idx));

% ----- 4. Časová smyčka IMEX Euler -----
I = eye(size(D2));
A = I - nu * dt * D2;         % implicitní difúze
N_steps = T / dt;

for step = 1:N_steps
    u_x = D1 * u;
    conv = u .* u_x;
    rhs = u - dt * conv;
    u = A \ rhs;
end

% ----- 5. Chyba -----
u_exact = @(x,t) (2*nu*pi*exp(-nu*pi^2*t).*sin(pi*x)) ./ (c + exp(-nu*pi^2*t).*cos(pi*x));
u_ref = u_exact(cheb_pts(idx), T);
err = u - u_ref;
max_err = norm(err, inf);
fprintf('Max. chyba v t = %.2f:  %.3e\n', T, max_err);

% ----- 6. Vykreslení -----
u_full = zeros(N, 1);
u_full(idx) = u;
u_cheb = chebfun(u_full, dom);
u_exact_cheb = chebfun(@(x) u_exact(x, T), dom);
err_cheb = u_cheb - u_exact_cheb;

figure;
subplot(1,2,1)
plot(u_cheb, 'b-', 'LineWidth', 1.5); hold on;
plot(u_exact_cheb, 'r--', 'LineWidth', 1.5);
xlabel('x'), ylabel('u(x,T)'), legend('IMEX Euler','Přesné'), title(sprintf('t = %.2f', T)), grid on

subplot(1,2,2)
plot(err_cheb, 'k-', 'LineWidth', 1.5);
xlabel('x'), ylabel('Chyba'), title(sprintf('Chyba (max = %.2e)', max_err)), grid on
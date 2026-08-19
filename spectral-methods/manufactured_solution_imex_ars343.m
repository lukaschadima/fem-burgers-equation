%% IMEX Runge-Kutta ARS(3,4,3) pro Burgersovu rovnici se zdrojovým členem
% Rovnice: u_t + u*u_x = nu * u_xx + f(x,t),   x∈[-1,1], t∈[0,T]
% Přesné řešení: u(x,t) = sin(πx)·sin(ωt)
% Samostatný skript – žádné externí soubory nejsou potřeba.

clear; close all;

% ============ PARAMETRY ============
nu      = 0.05;          % viskozita
omega   = 2*pi;          % frekvence zdroje
dom     = [-1, 1];       % výpočetní oblast
Tfinal  = 10.67;          % konečný čas
N       = 32;            % počet Čebyševových bodů
dt      = 1e-4;          % časový krok
tic;
% ============ ČEBYŠEVOVY BODY A DERIVAČNÍ MATICE ============
x   = chebpts(N, dom);
D1  = diffmat(N, 1, dom);       % 1. derivace
D2  = diffmat(N, 2, dom);       % 2. derivace

% ============ POČÁTEČNÍ PODMÍNKA (nulová) ============
U = zeros(N, 1);
U(1) = 0;   U(end) = 0;         % Dirichletovy okraje

% ============ ZDROJOVÝ ČLEN A PŘESNÉ ŘEŠENÍ ============
source_func = @(x,t) omega .* sin(pi*x) .* cos(omega*t) ...
               + (pi/2) .* sin(2*pi*x) .* sin(omega*t).^2 ...
               + nu*pi^2 .* sin(pi*x) .* sin(omega*t);

u_exact = @(x,t) sin(pi*x) .* sin(omega*t);

% ============ BUTCHEROVA TABULKA ARS(3,4,3) ============
gam = 0.4358665215484241;
b1  = -1.5*gam^2 + 4*gam - 0.25;
b2  =  1.5*gam^2 - 5*gam + 1.25;

Ai = [0,           0,          0,    0;
      0,           gam,        0,    0;
      0,  (1-gam)/2,            gam,  0;
      0,           b1,         b2,   gam];
bi = [0; b1; b2; gam];
ci = [0; gam; (1+gam)/2; 1];

a31 = 0.3212788860;   a32 = 0.3966543747;
a41 = -0.105858296;   a42 = 0.5529291479;   a43 = 0.5529291479;

Ae = [0,    0,    0,    0;
      gam,  0,    0,    0;
      a31,  a32,  0,    0;
      a41,  a42,  a43,  0];
be = [0; b1; b2; gam];
ce = [0; gam; (1+gam)/2; 1];

s = size(Ai,1);   % počet stupňů (4)

% ============ MATICE M (jednotková) a K (difúze) ============
I_mat = eye(N);
K = -D2;                     % tak, aby -ν*K = ν*D2

% ============ ČASOVÁ SMYČKA ============
t = 0;
nsteps = floor(Tfinal / dt);
fprintf('Spouštím ARS(3,4,3) s Δt = %.1e, počet kroků = %d\n', dt, nsteps);

for step = 1:nsteps
    U = step_imex_rk(U, t, dt, I_mat, K, D1, nu, source_func, ...
                     Ai, Ae, bi, be, ci, x);
    t = t + dt;
end
toc;
% ============ CHYBA VŮČI PŘESNÉMU ŘEŠENÍ ============
u_ref = u_exact(x, Tfinal);
err   = U - u_ref;
max_err = norm(err, inf);
fprintf('Max. chyba v t = %.2f:  %.3e\n', Tfinal, max_err);

% ============ VIZUALIZACE ============
u_cheb     = chebfun(U, dom);
u_ex_cheb  = chebfun(@(x) u_exact(x, Tfinal), dom);
err_cheb   = u_cheb - u_ex_cheb;

figure;
subplot(1,2,1)
plot(u_cheb, 'b-', 'LineWidth', 1.5); hold on
plot(u_ex_cheb, 'r--', 'LineWidth', 1.5)
xlabel('x'), ylabel('u(x,T)')
legend('IMEX RK ARS(3,4,3)', 'Přesné', 'Location', 'best')
title(sprintf('t = %.2f (Δt = %.1e)', Tfinal, dt))
grid on

subplot(1,2,2)
plot(err_cheb, 'k-', 'LineWidth', 1.5)
xlabel('x'), ylabel('Chyba')
title(sprintf('Chyba (max = %.2e)', max_err))
grid on

% ========================================================================
%  LOKÁLNÍ FUNKCE: JEDEN KROK IMEX RK
% ========================================================================
function U_new = step_imex_rk(U_prev, t_old, dt, M, K, D1, nu, source, ...
                               Ai, Ae, bi, be, ci, x)
% U_prev   ... řešení v čase t_old (sloupcový vektor délky N)
% t_old    ... počáteční čas kroku
% dt       ... časový krok
% M, K     ... hmotnostní a difúzní matice (N×N)
% D1       ... matice 1. derivace
% nu       ... viskozita
% source   ... zdrojová funkce @(x,t)
% Ai, Ae, bi, be, ci  ... Butcherovy tabulky (padded)
% x        ... Čebyševovy body (sloupcový vektor)

N_dof = length(U_prev);
s = size(Ai,1);               % počet stupňů

% Úložiště mezivýsledků
U_stages = zeros(N_dof, s);   % řešení na stupních
N_stages = zeros(N_dof, s);   % konvektivní členy N(U_j)
F_stages = zeros(N_dof, s);   % zdrojové členy F(t_j)

% ----- První stupeň (triviální, explicitní) -----
U_stages(:,1) = U_prev;
N_stages(:,1) = U_prev .* (D1 * U_prev);
F_stages(:,1) = source(x, t_old + ci(1)*dt);

% ----- Stupně 2 ... s -----
for i = 2:s
    t_i   = t_old + ci(i)*dt;
    Fi    = source(x, t_i);            % zdroj v čase stupně
    F_stages(:,i) = Fi;

    % Sestavení pravé strany
    rhs = M * U_prev;
    for j = 1:i-1
        % explicitní konvekce
        rhs = rhs + dt * Ae(i,j) * (-N_stages(:,j));
        % implicitní difuze + zdroj z předchozích stupňů
        rhs = rhs + dt * Ai(i,j) * (-nu * K * U_stages(:,j) + F_stages(:,j));
    end
    % implicitní zdroj aktuálního stupně
    rhs = rhs + dt * Ai(i,i) * Fi;

    % Matice soustavy
    A_mat = M + dt * Ai(i,i) * nu * K;

    % Dirichletovy okrajové podmínky (u = 0 na obou koncích)
    A_mat(1,:)   = 0;   A_mat(1,1)     = 1;   rhs(1)   = 0;
    A_mat(end,:) = 0;   A_mat(end,end) = 1;   rhs(end) = 0;

    U_stages(:,i) = A_mat \ rhs;
    N_stages(:,i) = U_stages(:,i) .* (D1 * U_stages(:,i));
end

% ----- Finální sestavení řešení v novém čase -----
rhs_final = M * U_prev;
for j = 1:s
    rhs_final = rhs_final ...
        + dt * bi(j) * (-nu * K * U_stages(:,j) + F_stages(:,j)) ...
        + dt * be(j) * (-N_stages(:,j));
end

A_final = M;    % hmotnostní matice (zde identita)
% Dirichlet pro finální soustavu
A_final(1,:)   = 0;   A_final(1,1)     = 1;   rhs_final(1)   = 0;
A_final(end,:) = 0;   A_final(end,end) = 1;   rhs_final(end) = 0;

U_new = A_final \ rhs_final;
end
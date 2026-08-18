function U = step_explicit_rk(p, n, x, U_prev, t_old, dt, params, T)
% U = step_explicit_rk(p, n, x, U_prev, t_old, dt, params, T)
%
% Jeden casovy krok explicitnim Runge-Kutta schematu.
% Aplikuje ERK na cely pravy clen: M*dU/dt = -N(U) - nu*K*U + F(t).
%
% POZNAMKA: Protoze difuze je zahrnuta explicitne, metoda je
% podmineně stabilni. Vyžaduje velmi male dt pro velke nu.
%
% Vstupy:
%   p, n, x, params  ... obvykle parametry (viz solve_burgers_equation)
%   T                ... struct z butcher_tables (typ 'explicit')
%
% Algoritmus:
%   Pro kazdy stupen i = 1,...,s:
%     U_i  = U_prev + dt * sum_{j<i} A(i,j) * K_j
%     M K_i = -N(U_i) - nu*K*U_i + F(t + c_i*dt)
%   U_new = U_prev + dt * sum_j b_j * K_j

nu     = params.nu;
s      = T.s;
N_dof  = length(U_prev);

% Hmotnostni matice M a tuhost K nezavisi na U ani t.
% Sestavime jednou z U_prev (jen abychom ziskali M, K).
[M, K_mat, ~, ~, ~] = assemble_burgers_system(p, n, x, U_prev, t_old, params);

% Modifikujeme M pro okrajove podminky: radek 1 a end -> identita
M_bc = M;
M_bc(1,:)   = 0;  M_bc(1,1)     = 1;
M_bc(end,:) = 0;  M_bc(end,end) = 1;

% LU rozklad M_bc jednou (M nezavisi na case ani reseni)
[L_M, U_M, P_M] = lu(M_bc);

% Ulozeni stupnovych derivaci K_i (jako vektory velikosti N_dof)
stage_deriv = zeros(N_dof, s);

for i = 1:s
    t_i = t_old + T.c(i) * dt;

    % --- Hodnota reseni na zacatku stupne ---
    U_i = U_prev;
    for j = 1:i-1
        U_i = U_i + dt * T.A(i,j) * stage_deriv(:,j);
    end
    % Okrajove podminky na U_i (pro jistotu)
    U_i(1) = 0;  U_i(end) = 0;

    % --- Pravy clen: -N(U_i) - nu*K*U_i + F(t_i) ---
    [~, K_loc, ~, N_i, F_i] = assemble_burgers_system(p, n, x, U_i, t_i, params);
    rhs_i = -N_i - nu * K_loc * U_i + F_i;

    % Okrajove podminky na pravy clen (derivace na hranici = 0)
    rhs_i(1)   = 0;
    rhs_i(end) = 0;

    % --- Reseni M * K_i = rhs_i (LU uz rozlozeno) ---
    stage_deriv(:,i) = U_M \ (L_M \ (P_M * rhs_i));
end

% --- Finalni update: U = U_prev + dt * sum_j b_j * K_j ---
U = U_prev;
for j = 1:s
    U = U + dt * T.b(j) * stage_deriv(:,j);
end
U(1) = 0;  U(end) = 0;
end

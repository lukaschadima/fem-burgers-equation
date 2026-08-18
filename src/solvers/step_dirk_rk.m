function U = step_dirk_rk(p, n, x, U_prev, t_old, dt, params, T)
% U = step_dirk_rk(p, n, x, U_prev, t_old, dt, params, T)
%
% Jeden casovy krok plne implicitnim SDIRK schematu.
% KONVEKCE I DIFUZE jsou implicitni --> Newtonova metoda na kazdem stupni.
%
% Na stupni i resime nelinearni soustavu:
%
%   R(U_i) = M*(U_i - rhs_explicit_i) / (dt*A(i,i))
%            + N(U_i) + nu*K*U_i - F(t_i) = 0
%
% ekvivalentne:
%
%   R(U_i) = M*U_i + dt*A(i,i)*(N(U_i) + nu*K*U_i - F(t_i)) - B_i = 0
%
% kde B_i = M*U_prev + dt*sum_{j<i} A(i,j)*(-N(U_j) - nu*K*U_j + F_j)
%
% Jacobian: J = M + dt*A(i,i)*(A_jac(U_i) + nu*K)
%
% Finalni update:
%   M*U_new = M*U_prev + dt*sum_j b(j)*(-N(U_j) - nu*K*U_j + F_j)
%
% Vstupy:
%   T ... struct z butcher_tables, T.type = 'dirk'

nu    = params.nu;
s     = T.s;
N_dof = length(U_prev);

% Sestaveni M a K (nezavisi na U ani t -- jen jednou)
[M, K_mat, ~, ~, ~] = assemble_burgers_system(p, n, x, U_prev, t_old, params);

% Uloziste hodnot reseni na stupnich a odpovidajicich pravych stran
U_stages = zeros(N_dof, s);
RHS_g    = zeros(N_dof, s);   % RHS_g(:,i) = -N(U_i) - nu*K*U_i + F_i

% =========================================================
%  STUPNOVA SMYCKA  (Newton na kazdem stupni)
% =========================================================
for i = 1:s
    t_i = t_old + T.c(i) * dt;

    % --- Explicitni cast prave strany z predchozich stupnu ---
    B_i = M * U_prev;
    for j = 1:i-1
        B_i = B_i + dt * T.A(i,j) * RHS_g(:,j);
    end

    % --- Pocatecni odhad pro Newton: linearizujeme difuzi ---
    % Pro prvni stupen: U_guess = U_prev
    % Pro dalsi: pouzijeme vysledek predchoziho stupne
    if i == 1
        U_i = U_prev;
    else
        U_i = U_stages(:,i-1);
    end
    U_i(1) = 0;  U_i(end) = 0;

    % --- Newtonova iterace ---
    for iter = 1:50
        [~, K_loc, A_jac, N_i, F_i] = assemble_burgers_system(p, n, x, U_i, t_i, params);

        % Reziduum: R = M*U_i + dt*aii*(N_i + nu*K*U_i - F_i) - B_i
        aii = T.A(i,i);
        R   = M*U_i + dt*aii*(N_i + nu*K_loc*U_i - F_i) - B_i;
        J   = M + dt*aii*(A_jac + nu*K_loc);

        % Okrajove podminky (Dirichlet: U=0 na hranici)
        R(1)   = U_i(1);    J(1,:)   = 0;  J(1,1)     = 1;
        R(end) = U_i(end);  J(end,:) = 0;  J(end,end) = 1;

        dU  = J \ (-R);
        U_i = U_i + dU;

        if norm(dU) < 1e-13, break; end
    end

    U_stages(:,i) = U_i;

    % Ulozeni -N(U_i) - nu*K*U_i + F_i pro finalni update
    [~, K_loc, ~, N_i, F_i] = assemble_burgers_system(p, n, x, U_i, t_i, params);
    RHS_g(:,i) = -N_i - nu*K_loc*U_i + F_i;
end

% =========================================================
%  FINALNI UPDATE
%  M*U_new = M*U_prev + dt * sum_j b(j) * RHS_g(:,j)
% =========================================================
rhs_final = M * U_prev;
for j = 1:s
    rhs_final = rhs_final + dt * T.b(j) * RHS_g(:,j);
end

% Okrajove podminky na M
M_bc = M;
M_bc(1,:)   = 0;  M_bc(1,1)     = 1;  rhs_final(1)   = 0;
M_bc(end,:) = 0;  M_bc(end,end) = 1;  rhs_final(end) = 0;

U = M_bc \ rhs_final;
end

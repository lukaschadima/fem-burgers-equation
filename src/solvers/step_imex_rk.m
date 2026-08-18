function U = step_imex_rk(p, n, x, U_prev, t_old, dt, params, T)
% U = step_imex_rk(p, n, x, U_prev, t_old, dt, params, T)
%
% Jeden časový krok IMEX Runge-Kutta schematu (padded ARS tabulka).
%
% Rozdělení pravé strany ODE M*dU/dt = f + g:
%   f(U) = -N(U)           ... explicitní část (konvekce, nelineární)
%   g(U) = -nu*K*U + F(t)  ... implicitní část (difuze + zdroj, lineární v U)
%
% Na každém stupni i se řeší LINEÁRNÍ soustava:
%   (M + dt*Aii(i,i)*nu*K) * U_i = rhs_i,
% kde rhs_i obsahuje:
%   - M*U_prev
%   - explicitní konvekci z předchozích stupňů (s koeficienty Ae)
%   - implicitní difuzi a zdroj z předchozích stupňů (s koeficienty Ai)
%   - zdrojový člen aktuálního stupně (s koeficientem Aii(i,i))
%
% Padded tabulka: první stupeň je triviální (Ai(1,:)=0, Ae(1,:)=0),
% proto U_1 = U_prev (řešíme jen M*U_1 = M*U_prev).
%
% Finalní update:
%   M*U_new = M*U_prev + dt * sum_j [ bi(j)*(-nu*K*U_j+F_j) + be(j)*(-N_j) ]

    nu    = params.nu;
    s     = T.s;           % počet stupňů (padded)
    N_dof = length(U_prev);

    % --- Sestavení M a K (nezávisí na U ani t) ---
    [M, K_mat, ~, ~, ~] = assemble_burgers_system(p, n, x, U_prev, t_old, params);

    % --- Úložiště hodnot na stupních ---
    U_stages = zeros(N_dof, s);   % U_stages(:,i) = řešení na i-tém stupni
    Nf       = zeros(N_dof, s);   % Nf(:,i)       = N(U_i)
    Ff       = zeros(N_dof, s);   % Ff(:,i)       = F(t_i)

    % =========================================================
    %  STUPŇOVÁ SMYČKA
    % =========================================================
    for i = 1:s
        t_i = t_old + T.ci(i) * dt;   % čas aktuálního stupně

        % --- Sestavení pravé strany stupně i ---
        rhs = M * U_prev;   % začínáme s historickým příspěvkem

        % Předchozí stupně (j < i)
        for j = 1:i-1
            % Explicitní část: Ae(i,j) * f(U_j)
            rhs = rhs + dt * T.Ae(i,j) * (-Nf(:,j));

            % Implicitní část: Ai(i,j) * g(U_j) = Ai(i,j) * (-nu*K*U_j + F_j)
            rhs = rhs + dt * T.Ai(i,j) * (-nu * K_mat * U_stages(:,j) + Ff(:,j));
        end

        % Zdrojový člen aktuálního stupně (pouze implicitní část, protože zdroj
        % byl zařazen do g; difuzní část -nu*K*U_i je na levé straně soustavy)
        [~, ~, ~, ~, F_i] = assemble_burgers_system(p, n, x, U_prev, t_i, params);
        Ff(:,i) = F_i;   % uložíme pro pozdější použití v dalších stupních a ve finále
        rhs = rhs + dt * T.Ai(i,i) * F_i;

        % --- Sestava levé strany: (M + dt*Ai(i,i)*nu*K) ---
        A_sys = M + dt * T.Ai(i,i) * nu * K_mat;

        % --- Okrajové podmínky (Dirichlet: U=0 na hranici) ---
        A_sys(1,:)   = 0;  A_sys(1,1)     = 1;  rhs(1)   = 0;
        A_sys(end,:) = 0;  A_sys(end,end) = 1;  rhs(end) = 0;

        % --- Řešení lineární soustavy (žádný Newton!) ---
        U_stages(:,i) = A_sys \ rhs;

        % --- Uložení nelineární konvekce pro další stupně ---
        [~, ~, ~, N_i, ~] = assemble_burgers_system(p, n, x, U_stages(:,i), t_i, params);
        Nf(:,i) = N_i;
    end

    % =========================================================
    %  FINÁLNÍ UPDATE
    %  M*U_new = M*U_prev + dt * sum_j [ bi(j)*(-nu*K*U_j+F_j) + be(j)*(-N_j) ]
    % =========================================================
    rhs_final = M * U_prev;
    for j = 1:s
        rhs_final = rhs_final ...
            + dt * T.bi(j) * (-nu * K_mat * U_stages(:,j) + Ff(:,j)) ...
            + dt * T.be(j) * (-Nf(:,j));
    end

    % Okrajové podmínky pro finální systém
    M_bc = M;
    M_bc(1,:)   = 0;  M_bc(1,1)     = 1;  rhs_final(1)   = 0;
    M_bc(end,:) = 0;  M_bc(end,end) = 1;  rhs_final(end) = 0;

    U = M_bc \ rhs_final;
end
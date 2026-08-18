function [x, U_history, time_steps] = solve_burgers_equation(p, n, dt, params, method)
    % [x, U_history, time_steps] = solve_burgers_equation(p, n, dt, params, method)
    %
    % method: 'implicit_euler' | 'explicit_euler' | 'crank_nicolson' | 'imex_euler'

    L         = params.L;
    num_nodes = n*p + 1;

    x = zeros(num_nodes, 1);
    for e = 1:n
        x_left = (e-1)*L/n;
        for i = 0:p
            x((e-1)*p + i + 1) = x_left + i*(L/n)/p;
        end
    end

    num_steps  = ceil(params.T_final / dt);
    time_steps = linspace(0, params.T_final, num_steps + 1);
    dt_eff     = time_steps(2) - time_steps(1);

    U_history      = zeros(num_nodes, num_steps + 1);
    U_history(:,1) = params.ic(x);

    for k = 1:num_steps
        t_old  = time_steps(k);
        t_new  = time_steps(k+1);
        U_prev = U_history(:,k);

        switch method
            % --- Puvodni jednoduche metody (zachovany beze zmeny) ---
            case 'implicit_euler'
                U_history(:,k+1) = step_implicit_euler(p, n, x, U_prev, t_new, dt_eff, params);
            case 'explicit_euler'
                U_history(:,k+1) = step_explicit_euler(p, n, x, U_prev, t_old, dt_eff, params);
            case 'crank_nicolson'
                U_history(:,k+1) = step_crank_nicolson(p, n, x, U_prev, t_old, t_new, dt_eff, params);
            case 'imex_euler'
                U_history(:,k+1) = step_imex_euler(p, n, x, U_prev, t_new, dt_eff, params);

            % --- Explicitni Runge-Kutta (erk1 = forward Euler) ---
            case {'erk1', 'erk2', 'erk3', 'erk4'}
                T = butcher_tables(method);
                U_history(:,k+1) = step_explicit_rk(p, n, x, U_prev, t_old, dt_eff, params, T);

            % --- Plne implicitni SDIRK (Newton na kazdem stupni) ---
            case {'dirk2', 'dirk3'}
                T_rk = butcher_tables(method);
                U_history(:,k+1) = step_dirk_rk(p, n, x, U_prev, t_old, dt_eff, params, T_rk);

            % --- IMEX-ARS Runge-Kutta schemata ---
            case {'imex_ars111', 'imex_ars122', 'imex_ars232', 'imex_ars343'}
                T = butcher_tables(method);
                U_history(:,k+1) = step_imex_rk(p, n, x, U_prev, t_old, dt_eff, params, T);

            otherwise
                error(['Neznamá metoda: "%s".\n' ...
                       'Dostupne: implicit_euler, explicit_euler, crank_nicolson, imex_euler,\n' ...
                       '          erk1, erk2, erk3, erk4, dirk2, dirk3,\n' ...
                       '          imex_ars111, imex_ars122, imex_ars232, imex_ars343'], method);
        end
    end
end

% =========================================================================

function U = step_implicit_euler(p, n, x, U_prev, t, dt, params)
    % R = M*(U-U_prev)/dt + N(U) + nu*K*U - F(t) = 0
    % J = M/dt + A_jac + nu*K
    nu = params.nu;
    U  = U_prev;

    for iter = 1:50
        [M, K, A_jac, N_vec, F_vec] = assemble_burgers_system(p, n, x, U, t, params);

        R = M*(U - U_prev)/dt + N_vec + nu*K*U - F_vec;
        J = M/dt + A_jac + nu*K;

        J(1,:)   = 0; J(1,1)     = 1; R(1)   = U(1);
        J(end,:) = 0; J(end,end) = 1; R(end) = U(end);

        dU = J \ (-R);
        U  = U + dU;
        if norm(dU) < 1e-14, break; end
    end
end

function U = step_explicit_euler(p, n, x, U_prev, t, dt, params)
    % M*U_new = M*U_prev + dt*(F(t) - N(U_prev) - nu*K*U_prev)
    nu = params.nu;
    [M, K, ~, N_vec, F_vec] = assemble_burgers_system(p, n, x, U_prev, t, params);

    rhs = M*U_prev + dt*(F_vec - N_vec - nu*K*U_prev);

    M(1,:)   = 0; M(1,1)     = 1; rhs(1)   = 0;
    M(end,:) = 0; M(end,end) = 1; rhs(end) = 0;

    U = M \ rhs;
end

function U = step_crank_nicolson(p, n, x, U_prev, t_old, t_new, dt, params)
    % R = M*(U-U_prev)/dt + 0.5*(N_new+N_prev) + 0.5*nu*K*(U+U_prev) - 0.5*(F_new+F_prev)
    % J = M/dt + 0.5*(A_jac + nu*K)
    nu = params.nu;
    [M, K, ~, N_prev, F_prev] = assemble_burgers_system(p, n, x, U_prev, t_old, params);
    U = U_prev;

    for iter = 1:50
        [~, ~, A_jac, N_new, F_new] = assemble_burgers_system(p, n, x, U, t_new, params);

        R = M*(U - U_prev)/dt ...
            + 0.5*(N_new + N_prev) ...
            + 0.5*nu*K*(U + U_prev) ...
            - 0.5*(F_new + F_prev);
        J = M/dt + 0.5*(A_jac + nu*K);

        J(1,:)   = 0; J(1,1)     = 1; R(1)   = U(1);
        J(end,:) = 0; J(end,end) = 1; R(end) = U(end);

        dU = J \ (-R);
        U  = U + dU;
        if norm(dU) < 1e-14, break; end
    end
end

function U = step_imex_euler(p, n, x, U_prev, t_new, dt, params)
    % (M/dt + nu*K)*U_new = M/dt*U_prev + F(t_new) - N(U_prev)
    % N explicitne (U_prev), difuze implicitne - linearni soustava
    nu = params.nu;
    [M, K, ~, N_vec, F_vec] = assemble_burgers_system(p, n, x, U_prev, t_new, params);

    A   = M/dt + nu*K;
    rhs = M/dt * U_prev + F_vec - N_vec;

    A(1,:)   = 0; A(1,1)     = 1; rhs(1)   = 0;
    A(end,:) = 0; A(end,end) = 1; rhs(end) = 0;

    U = A \ rhs;
end

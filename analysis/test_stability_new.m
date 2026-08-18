function test_stability_new()
% Srovnani stability a casove konvergence vsech casovych metod.
% Example 2 (Manufactured, sin*sin) -- zname presne reseni.
%
% Sada metod:
%   Explicitni RK     : explicit_euler, erk2, erk3, erk4
%   Plne implicitni   : implicit_euler, dirk2, dirk3
%   IMEX              : imex_euler (=ARS111), imex_ars122, imex_ars232, imex_ars343
%   Dalsi             : crank_nicolson
%
% Pro male dt vsechny metody konverguji -- vidime rad.
% Pro velke dt explicitni RK diverguji -- ukazujeme nestabilitu.

    p      = 3;
    n      = 64;
    nu     = 0.05;
    omega  = 1;
    params = burgers_params(2, nu, omega, 1.0);

    % CFL odhad pro FEM (konzervativnejsi nez FD kvuli konzistentni M)
    h        = params.L / n;
    dt_crit  = h^2 / (6*nu);   % teoreticka hranice z FEM analyzy
    dt_fd    = h^2 / (2*nu);   % FD analogon (pro srovnani)

    fprintf('=== STABILITA + CASOVA KONVERGENCE ===\n');
    fprintf('p=%d, n=%d, h=%.4f, nu=%.3f, omega=%d, T=%.2f\n', p, n, h, nu, omega, params.T_final);
    fprintf('Teoreticky dt_crit (FEM) = %.5f\n', dt_crit);
    fprintf('Analogon FD             = %.5f\n\n', dt_fd);

    % dt hodnoty: od velkych (nestabilita ERK) po male (konvergence)
    dt_values = [ 0.02, 0.01, 0.005, 0.002, 0.001];

    % =========================================================
    %  SKUPINY METOD  (pro prehlednost tabulky a grafu)
    %  Kazda skupina = cell array {method_string, label, color, marker}
    % =========================================================
    method_defs = {
        % -- Explicitni RK (budou nestabilni pro velke dt) --
        %'explicit_euler',  'ERK1 (Expl.Euler)', [0.8 0.2 0.2], 'v';
        %'erk2',            'ERK2 (Heun)',        [0.9 0.4 0.0], 'd';
        % 'erk3',            'ERK3',               [0.8 0.6 0.0], 'p';
        % 'erk4',            'ERK4',               [0.6 0.4 0.0], 'h';
        % -- Plne implicitni (bezpodm. stabilni, drahe) --
        'implicit_euler',  'DIRK1 (Impl.Euler)', [0.0 0.4 0.8], 'o';
        'dirk2',           'DIRK2',              [0.0 0.2 0.6], 's';
        'dirk3',           'DIRK3',              [0.0 0.0 0.4], '^';
        % -- IMEX (bezpodm. stabilni na difuzi, levne) --
        'imex_euler',      'IMEX-ARS(1,1,1)',    [0.0 0.6 0.3], 'o';
        'imex_ars122',     'IMEX-ARS(1,2,2)',    [0.0 0.5 0.0], 's';
        'imex_ars232',     'IMEX-ARS(2,3,2)',    [0.2 0.7 0.2], '^';
        'imex_ars343',     'IMEX-ARS(3,4,3)',    [0.4 0.8 0.4], 'd';
        % -- Ostatni --
        'crank_nicolson',  'Crank-Nicolson',     [0.5 0.0 0.5], '*';
    };

    n_methods = size(method_defs, 1);
    n_dt      = length(dt_values);

    fmt_hdr = '%-22s  %-8s  %-12s  %-12s  %-10s  %s\n';
    fmt_row = '%-22s  %-8.4f  %-12.3e  %-12.3e  %-10.3f  %s\n';
    fmt_div = '%-22s  %-8.4f  %-12s  %-12s  %-10s  NE\n';
    fmt_sep = '%-22s  %-8s  %-12s  %-12s  %-10s  %s\n';

    fprintf(fmt_hdr, 'Metoda', 'dt', 'L2 chyba', 'H1 chyba', 'Cas [s]', 'Stabilni');
    fprintf('%s\n', repmat('-', 1, 82));

    all_results = cell(n_methods, 1);

    for mi = 1:n_methods
        method = method_defs{mi, 1};
        label  = method_defs{mi, 2};

        res = struct('dt',     dt_values, ...
                     'L2',     nan(size(dt_values)), ...
                     'H1',     nan(size(dt_values)), ...
                     't_cpu',  nan(size(dt_values)), ...
                     'stable', false(size(dt_values)));

        for di = 1:n_dt
            dt = dt_values(di);
            try
                t_start = tic;
                [x, U_history, time_steps] = ...
                    solve_burgers_equation(p, n, dt, params, method);
                t_cpu = toc(t_start);

                U_final = U_history(:, end);
                stable  = isfinite(norm(U_final)) && norm(U_final) < 1e6;

                if stable
                    [L2, H1] = compute_errors_local(p, x, U_final, time_steps(end), params);
                    res.L2(di)     = L2;
                    res.H1(di)     = H1;
                    res.t_cpu(di)  = t_cpu;
                    res.stable(di) = true;
                    fprintf(fmt_row, label, dt, L2, H1, t_cpu, 'ANO');
                else
                    res.t_cpu(di) = t_cpu;
                    fprintf(fmt_div, label, dt, 'Inf', 'Inf', sprintf('%.3f', t_cpu));
                end
            catch ME
                fprintf(fmt_div, label, dt, 'ERR', 'ERR', '--');
            end
        end

        % Casovy rad konvergence (jen ze stabilnich bodu)
        stable_idx = find(res.stable);
        if length(stable_idx) >= 2
            dts  = dt_values(stable_idx);
            errs = res.L2(stable_idx);
            ord  = mean(diff(log(errs)) ./ diff(log(dts)));
            fprintf('  -> prumerny casovy rad L2 (%s): %.2f\n', label, ord);
        end
        fprintf('\n');

        all_results{mi} = res;
    end

    % =========================================================
    %  GRAFY
    % =========================================================

    % --- Graf 1: L2 chyba vs dt (casova konvergence) ---
    fig1 = figure('Position', [50 50 720 500]);
    hold on; grid on;
    for mi = 1:n_methods
        res    = all_results{mi};
        label  = method_defs{mi, 2};
        col    = method_defs{mi, 3};
        mark   = method_defs{mi, 4};
        mask   = res.stable;
        if sum(mask) >= 2
            loglog(res.dt(mask), res.L2(mask), ...
                'Color', col, 'Marker', mark, 'LineWidth', 1.8, ...
                'MarkerSize', 7, 'DisplayName', label);
        end
    end
    % Referencni sklony
    dt_ref = dt_values(end-2:end);
    loglog(dt_ref, 2e-3*(dt_ref/dt_ref(1)).^1, 'k--', 'LineWidth', 1, ...
           'DisplayName', 'O(\Deltat)');
    loglog(dt_ref, 8e-5*(dt_ref/dt_ref(1)).^2, 'k-.',  'LineWidth', 1, ...
           'DisplayName', 'O(\Deltat^2)');
    loglog(dt_ref, 3e-6*(dt_ref/dt_ref(1)).^3, 'k:',   'LineWidth', 1, ...
           'DisplayName', 'O(\Deltat^3)');
    set(gca, 'XScale', 'log', 'YScale', 'log');
    xlabel('\Deltat', 'FontSize', 12);
    ylabel('||u - u_h||_{L^2}', 'FontSize', 12);
    title(sprintf('Casova konvergence L^2  (p=%d, n=%d, \\nu=%.3f)', p, n, nu));
    legend('Location', 'southeast', 'FontSize', 8);
    exportgraphics(fig1, 'stability_L2_all.pdf', 'ContentType', 'vector');

    % --- Graf 2: CPU cas vs dt ---
    fig2 = figure('Position', [50 50 720 500]);
    hold on; grid on;
    for mi = 1:n_methods
        res   = all_results{mi};
        label = method_defs{mi, 2};
        col   = method_defs{mi, 3};
        mark  = method_defs{mi, 4};
        mask  = res.stable & ~isnan(res.t_cpu);
        if any(mask)
            loglog(res.dt(mask), res.t_cpu(mask), ...
                'Color', col, 'Marker', mark, 'LineWidth', 1.8, ...
                'MarkerSize', 7, 'DisplayName', label);
        end
    end
    set(gca, 'XScale', 'log', 'YScale', 'log');
    xlabel('\Deltat', 'FontSize', 12);
    ylabel('CPU cas [s]', 'FontSize', 12);
    title(sprintf('Vypocetni cas  (p=%d, n=%d, \\nu=%.3f)', p, n, nu));
    legend('Location', 'best', 'FontSize', 8);
    exportgraphics(fig2, 'stability_cpu_all.pdf', 'ContentType', 'vector');

    % --- Graf 3: Efektivita (L2 vs CPU) -- Pareto porovnani ---
    fig3 = figure('Position', [50 50 720 500]);
    hold on; grid on;
    for mi = 1:n_methods
        res   = all_results{mi};
        label = method_defs{mi, 2};
        col   = method_defs{mi, 3};
        mark  = method_defs{mi, 4};
        mask  = res.stable & ~isnan(res.t_cpu) & ~isnan(res.L2);
        if sum(mask) >= 2
            loglog(res.t_cpu(mask), res.L2(mask), ...
                'Color', col, 'Marker', mark, 'LineWidth', 1.8, ...
                'MarkerSize', 7, 'DisplayName', label);
        end
    end
    set(gca, 'XScale', 'log', 'YScale', 'log');
    xlabel('CPU cas [s]', 'FontSize', 12);
    ylabel('||u - u_h||_{L^2}', 'FontSize', 12);
    title('Efektivita: chyba vs vypocetni cas (mene-vlevo-dole = lepsi)');
    legend('Location', 'northeast', 'FontSize', 8);
    exportgraphics(fig3, 'stability_efficiency.pdf', 'ContentType', 'vector');
end

% -------------------------------------------------------------------------

function [L2_err, H1_err] = compute_errors_local(p, x, U, t, params)
    n_elem = (length(x) - 1) / p;
    [gauss_pts, gauss_wts] = gauss_quadrature_info(max(3, p+2));
    L2_sq = 0; H1_sq = 0;

    for e = 1:n_elem
        s      = (e-1)*p + 1;
        x_l    = x(s); x_r = x(s+p);
        u_elem = U(s:s+p);
        h      = x_r - x_l; Jac = h/2;

        for k = 1:length(gauss_pts)
            xi = gauss_pts(k); w = gauss_wts(k);
            [phi, dphi_dxi] = lagrange_basis(p, xi);
            x_phys = (x_l+x_r)/2 + Jac*xi;
            [u_ex, ux_ex] = params.exact(x_phys, t);
            u_fem  = dot(u_elem, phi);
            ux_fem = dot(u_elem, dphi_dxi) * (2/h);
            L2_sq  = L2_sq + w*(u_ex - u_fem)^2  * Jac;
            H1_sq  = H1_sq + w*(ux_ex - ux_fem)^2 * Jac;
        end
    end

    L2_err = sqrt(L2_sq);
    H1_err = sqrt(H1_sq);
end

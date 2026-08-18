function test_stability_legacy()
    % Srovnani stability a vypocetniho casu
    % implicit_euler vs imex_euler vs explicit_euler
    % Example 2 (Manufactured, sin*sin)

    p      = 2;
    n      = 32;
    nu     = 0.05;
    omega  = 10;
    params = burgers_params(2, nu, omega, 1.0);

    % Teoreticky odhad dt_crit pro FEM je konzervativnejsi nez pro FD
    % Skutecna hranice lezi mezi 0.0005 a 0.001
    dt_crit = (params.L/n)^2 / (20*nu);

    dt_values = [ 0.005, 0.002, 0.001,0.0005, 0.0002, 0.0001];
    methods   = {'implicit_euler', 'imex_euler', 'explicit_euler'};
    labels    = {'Implicit Euler', 'IMEX Euler', 'Explicit Euler'};

    fprintf('=== STABILITA A CAS (p=%d, n=%d, nu=%.3f, omega=%d) ===\n', p, n, nu, omega);
    fprintf('h = %.4f,  dt_crit ≈ %.4f\n\n', params.L/n, dt_crit);

    fmt_hdr = '%-16s  %-8s  %-12s  %-12s  %-10s  %s\n';
    fmt_row = '%-16s  %-8.4f  %-12.3e  %-12.3e  %-10.3f  %s\n';
    fmt_div = '%-16s  %-8.4f  %-12s  %-12s  %-10s  NE\n';
    fprintf(fmt_hdr, 'Metoda', 'dt', 'L2 chyba', 'H1 chyba', 'Cas [s]', 'Stabilni');
    fprintf('%s\n', repmat('-', 1, 76));

    all_results = cell(length(methods), 1);
    colors      = {'b', 'g', 'r'};
    markers     = {'o', 's', '^'};

    for mi = 1:length(methods)
        method = methods{mi};
        res    = struct('dt', dt_values, 'L2', nan(size(dt_values)), ...
                        'H1', nan(size(dt_values)), 't_cpu', nan(size(dt_values)), ...
                        'stable', false(size(dt_values)));

        for di = 1:length(dt_values)
            dt = dt_values(di);
            try
                t_start = tic;
                [x, U_history, time_steps] = ...
                    solve_burgers_equation(p, n, dt, params, method);
                t_cpu = toc(t_start);

                U_final = U_history(:,end);
                stable  = isfinite(norm(U_final)) && norm(U_final) < 1e6;

                if stable
                    [L2, H1] = compute_errors_local(p, x, U_final, time_steps(end), params);
                    res.L2(di) = L2; res.H1(di) = H1;
                    res.t_cpu(di) = t_cpu; res.stable(di) = true;
                    fprintf(fmt_row, labels{mi}, dt, L2, H1, t_cpu, 'ANO');
                else
                    res.t_cpu(di) = t_cpu;
                    fprintf(fmt_div, labels{mi}, dt, 'Inf', 'Inf', sprintf('%.3f', t_cpu));
                end
            catch
                fprintf(fmt_div, labels{mi}, dt, 'Inf', 'Inf', '--');
            end
        end

        if sum(res.stable) >= 2
            dts  = dt_values(res.stable);
            errs = res.L2(res.stable);
            ord  = mean(diff(log(errs)) ./ diff(log(dts)));
            fprintf('  -> prumerny casovy rad L2: %.2f\n', ord);
        end

        all_results{mi} = res;
        fprintf('\n');
    end

    ttl = sprintf('\\nu=%.3f, \\omega=%d, n=%d, p=%d', nu, omega, n, p);

    % Graf 1: L2 chyba vs dt
    fig1 = figure('Position', [50 50 640 480]);
    hold on; grid on;
    for mi = 1:length(methods)
        res  = all_results{mi};
        mask = res.stable;
        if any(mask)
            loglog(res.dt(mask), res.L2(mask), ...
                [markers{mi} colors{mi} '-'], 'LineWidth', 2, ...
                'MarkerSize', 8, 'DisplayName', labels{mi});
        end
    end
    dt_ref = dt_values(end-2:end);
    %loglog(dt_ref, 0.05*dt_ref/dt_ref(1),          'k--', 'DisplayName', 'O(\Deltat)');
    %loglog(dt_ref, 0.001*(dt_ref/dt_ref(1)).^2,     'k:',  'DisplayName', 'O(\Deltat^2)');
    set(gca, 'XScale', 'log', 'YScale', 'log');
    xlabel('\Deltat'); ylabel('||u - u_h||_{L^2}');
    title(['Casova konvergence L^2']);
    legend('Location', 'southeast', 'FontSize', 9); grid on;
    exportgraphics(fig1, 'stability_L2.pdf', 'ContentType', 'vector');

    % Graf 2: CPU cas vs dt (obe osy log)
    fig2 = figure('Position', [50 50 640 480]);
    hold on; grid on;
    for mi = 1:length(methods)
        res  = all_results{mi};
        mask = res.stable & ~isnan(res.t_cpu);
        if any(mask)
            loglog(res.dt(mask), res.t_cpu(mask), ...
                [markers{mi} colors{mi} '-'], 'LineWidth', 2, ...
                'MarkerSize', 8, 'DisplayName', labels{mi});
        end
    end
    set(gca, 'XScale', 'log', 'YScale', 'log');
    xlabel('\Deltat'); ylabel('CPU cas [s]');
    title(['Výpocetní cas']);
    legend('Location', 'best', 'FontSize', 9); grid on;
    exportgraphics(fig2, 'stability_cpu.pdf', 'ContentType', 'vector');
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

function test_spatial_convergence()
    % Prostorova konvergence - implicit Euler, Example 1 (Cole-Hopf)
    % Maly dt => casova chyba zanedbatelna, vidime ciste prostorove rady

    nu       = 0.1;
    dt       = 5e-4;
    p_values = [1, 2, 3];
    n_values = [4, 8, 16, 32];
    params   = burgers_params(1, nu, 0.5);

    L2_errors = zeros(length(n_values), length(p_values));
    H1_errors = zeros(length(n_values), length(p_values));

    fprintf('=== PROSTOROVA KONVERGENCE (Implicit Euler, Cole-Hopf) ===\n');
    fprintf('nu=%.3f, T=%.2f, dt=%.0e\n\n', nu, params.T_final, dt);
    fprintf('%-5s  %-5s  %-8s  %-12s  %-6s  %-12s  %-6s\n', ...
        'p', 'n', 'h', 'L2 chyba', 'rad', 'H1 chyba', 'rad');
    fprintf('%s\n', repmat('-', 1, 62));

    for pi = 1:length(p_values)
        p = p_values(pi);
        for ni = 1:length(n_values)
            n = n_values(ni);
            h = params.L / n;

            [x, U_history, time_steps] = ...
                solve_burgers_equation(p, n, dt, params, 'crank_nicolson');

            [L2, H1] = compute_errors_local(p, x, U_history(:,end), time_steps(end), params);
            L2_errors(ni, pi) = L2;
            H1_errors(ni, pi) = H1;

            if ni == 1
                fprintf('%-5d  %-5d  %-8.4f  %-12.3e  %-6s  %-12.3e  %-6s\n', ...
                    p, n, h, L2, '--', H1, '--');
            else
                ord_L2 = log(L2_errors(ni-1,pi)/L2) / log(2);
                ord_H1 = log(H1_errors(ni-1,pi)/H1) / log(2);
                fprintf('%-5d  %-5d  %-8.4f  %-12.3e  %-6.2f  %-12.3e  %-6.2f\n', ...
                    p, n, h, L2, ord_L2, H1, ord_H1);
            end
        end
        fprintf('\n');
    end

    h_values = params.L ./ n_values;
    colors   = {'b', 'r', 'g'};
    markers  = {'o', 's', '^'};
    ttl      = sprintf('Cole-Hopf, \\nu=%.1f', nu);

    % Graf L2
    fig1 = figure('Position', [50 50 640 480]);
    hold on; grid on;
    for pi = 1:length(p_values)
        p   = p_values(pi);
        ref = L2_errors(1,pi) * (h_values/h_values(1)).^(p+1);
        loglog(h_values, L2_errors(:,pi), [markers{pi} colors{pi} '-'], ...
            'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', sprintf('p=%d', p));
        loglog(h_values, ref, [colors{pi} ':'], 'LineWidth', 1.2, 'HandleVisibility', 'off');
    end
    h_sorted = fliplr(h_values);
    h_labels = arrayfun(@(h) sprintf('1/%d',round(1/h)), h_sorted, 'UniformOutput', false);
    set(gca, 'XScale', 'log', 'YScale', 'log', ...
        'XTick', h_sorted, 'XTickLabel', h_labels);
    xlabel('h'); ylabel('||u - u_h||_{L^2}');
    title(['Prostorová konvergence L^2']);
    legend('Location', 'southeast'); grid on;
    exportgraphics(fig1, 'convergence_L2.pdf', 'ContentType', 'vector');

    % Graf H1
    fig2 = figure('Position', [50 50 640 480]);
    hold on; grid on;
    for pi = 1:length(p_values)
        p   = p_values(pi);
        ref = H1_errors(1,pi) * (h_values/h_values(1)).^p;
        loglog(h_values, H1_errors(:,pi), [markers{pi} colors{pi} '-'], ...
            'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', sprintf('p=%d', p));
        loglog(h_values, ref, [colors{pi} ':'], 'LineWidth', 1.2, 'HandleVisibility', 'off');
    end
    set(gca, 'XScale', 'log', 'YScale', 'log', ...
        'XTick', h_sorted, 'XTickLabel', h_labels);
    xlabel('h'); ylabel('|u - u_h|_{H^1}');
    title(['Prostorová konvergence H^1']);
    legend('Location', 'southeast'); grid on;
    exportgraphics(fig2, 'convergence_H1.pdf', 'ContentType', 'vector');
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

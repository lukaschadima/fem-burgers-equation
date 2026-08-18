function results = burgers_analysis(params, p, n, dt_values, methods)
    % results = burgers_analysis(params, p, n, dt_values, methods)
    %
    % Spusti vsechny kombinace dt x method, vypise tabulku chyb a stability.
    %
    % Priklad:
    %   params  = burgers_params(2, 0.05, 10);
    %   results = burgers_analysis(params, 2, 32, [0.05 0.01 0.005], ...
    %               {'implicit_euler','explicit_euler','crank_nicolson','imex_euler'});

    if nargin < 5
        methods = {'implicit_euler', 'explicit_euler', 'crank_nicolson', 'imex_euler'};
    end
    if ischar(methods)
        methods = {methods};
    end

    fprintf('=== Burgersova rovnice - analyza ===\n');
    fprintf('Example %d, p=%d, n=%d, nu=%.4f, T=%.3f\n\n', ...
        params.example_id, p, n, params.nu, params.T_final);

    fmt_hdr = '%-20s  %-10s  %-12s  %-12s  %s\n';
    fmt_row = '%-20s  %-10.5f  %-12.3e  %-12.3e  %s\n';
    fmt_div = '%-20s  %-10.5f  %-12s  %-12s  NE\n';

    fprintf(fmt_hdr, 'Metoda', 'dt', 'L2 chyba', 'H1 chyba', 'Stabilni');
    fprintf('%s\n', repmat('-', 1, 68));

    n_methods = length(methods);
    n_dt      = length(dt_values);
    results   = struct('method',{}, 'dt',{}, 'L2_error',{}, 'H1_error',{}, 'stable',{});
    idx       = 0;

    for m = 1:n_methods
        for d = 1:n_dt
            idx = idx + 1;
            dt  = dt_values(d);

            try
                [x, U_history, time_steps] = ...
                    solve_burgers_equation(p, n, dt, params, methods{m});

                U_final = U_history(:, end);
                stable  = isfinite(norm(U_final)) && norm(U_final) < 1e6;

                if stable
                    [L2_err, H1_err] = compute_errors(p, x, U_final, time_steps(end), params);
                else
                    L2_err = Inf; H1_err = Inf;
                end
            catch
                stable = false; L2_err = Inf; H1_err = Inf;
            end

            results(idx).method   = methods{m};
            results(idx).dt       = dt;
            results(idx).L2_error = L2_err;
            results(idx).H1_error = H1_err;
            results(idx).stable   = stable;

            if stable
                fprintf(fmt_row, methods{m}, dt, L2_err, H1_err, 'ANO');
            else
                fprintf(fmt_div, methods{m}, dt, 'Inf','Inf');
            end
        end

        if n_dt > 1
            % Konvergencni rad pro tuto metodu (jen stabilni body)
            stable_mask = [results(idx-n_dt+1:idx).stable];
            if sum(stable_mask) >= 2
                dts  = dt_values(stable_mask);
                errs = [results(idx-n_dt+1:idx).L2_error];
                errs = errs(stable_mask);
                ord  = mean(diff(log(errs)) ./ diff(log(dts)));
                fprintf('  -> prumerny casovy rad L2: %.2f\n', ord);
            end
        end

        fprintf('\n');
    end
end

% -------------------------------------------------------------------------

function [L2_err, H1_err] = compute_errors(p, x, U, t, params)
    n_elem    = (length(x) - 1) / p;
    quad_ord  = max(3, p+2);
    [gauss_pts, gauss_wts] = gauss_quadrature_info(quad_ord);

    L2_sq = 0;
    H1_sq = 0;

    for e = 1:n_elem
        s       = (e-1)*p + 1;
        x_left  = x(s);
        x_right = x(s+p);
        u_elem  = U(s:s+p);
        h       = x_right - x_left;
        Jac     = h/2;

        for k = 1:length(gauss_pts)
            xi = gauss_pts(k);
            w  = gauss_wts(k);

            [phi, dphi_dxi] = lagrange_basis(p, xi);
            x_phys = (x_left + x_right)/2 + Jac*xi;

            [u_ex, ux_ex] = params.exact(x_phys, t);
            u_fem  = dot(u_elem, phi);
            ux_fem = dot(u_elem, dphi_dxi) * (2/h);

            L2_sq = L2_sq + w*(u_ex - u_fem)^2  * Jac;
            H1_sq = H1_sq + w*(ux_ex - ux_fem)^2 * Jac;
        end
    end

    L2_err = sqrt(L2_sq);
    H1_err = sqrt(H1_sq);
end

function params = burgers_params(example_id, nu, varargin)
    % burgers_params(1, nu)              - Cole-Hopf, T=1.0
    % burgers_params(1, nu, T_final)     - Cole-Hopf, vlastni T
    % burgers_params(2, nu, omega)       - Manufactured, T=2pi/omega
    % burgers_params(2, nu, omega, T_final)

    params.example_id = example_id;
    params.L  = 1.0;
    params.nu = nu;
    params.c  = 2.0;  % Cole-Hopf konstanta, c > 1 zajistuje phi > 0

    switch example_id
        case 1
            params.omega   = 0;
            params.T_final = 1.0;
            if ~isempty(varargin)
                params.T_final = varargin{1};
            end
        case 2
            if isempty(varargin)
                error('Example 2 vyzaduje omega jako treti argument');
            end
            params.omega   = varargin{1};
            params.T_final = 2*pi / params.omega;
            if length(varargin) >= 2
                params.T_final = varargin{2};
            end
        otherwise
            error('example_id musi byt 1 nebo 2');
    end

    p_copy        = params;
    params.ic     = @(x)    ic_fun(x, p_copy);
    params.source = @(x, t) source_fun(x, t, p_copy);
    params.exact  = @(x, t) exact_fun(x, t, p_copy);
end

% -------------------------------------------------------------------------

function u0 = ic_fun(x, p)
    switch p.example_id
        case 1  % u(x,0) = 2*nu*pi*sin(pi*x) / (2 + cos(pi*x))
            u0 = 2*p.nu*pi * sin(pi*x) ./ (p.c + cos(pi*x));
        case 2  % u(x,0) = 0
            u0 = zeros(size(x));
    end
end

function f = source_fun(x, t, p)
    switch p.example_id
        case 1
            f = zeros(size(x));
        case 2
            % u = sin(pi*x)*sin(omega*t)
            % f = u_t + u*u_x - nu*u_xx
            w = p.omega;
            f = w*sin(pi*x).*cos(w*t) ...
                + (pi/2)*sin(2*pi*x).*sin(w*t).^2 ...
                + p.nu*pi^2*sin(pi*x).*sin(w*t);
    end
end

function [u, ux] = exact_fun(x, t, p)
    switch p.example_id
        case 1
            % u = 2*nu*pi * e*sin(pi*x) / (c + e*cos(pi*x))
            % kde e = exp(-nu*pi^2*t)
            e     = exp(-p.nu*pi^2*t);
            denom = p.c + e*cos(pi*x);
            u     = 2*p.nu*pi * e*sin(pi*x)             ./ denom;
            ux    = 2*p.nu*pi^2 * e * (p.c*cos(pi*x)+e) ./ denom.^2;
        case 2
            w  = p.omega;
            u  = sin(pi*x)    .* sin(w*t);
            ux = pi*cos(pi*x) .* sin(w*t);
    end
end

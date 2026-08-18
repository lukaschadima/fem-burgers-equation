function [M, K, A_jac, N_vec, F_vec] = assemble_burgers_system(p, n, x, U, t, params)
    % Vraci:
    %   M     - hmotnostni matice
    %   K     - difuzni matice (BEZ faktoru nu)
    %   A_jac - Jacobian nelinearniho clenu dN/dU
    %   N_vec - vektor nelinearniho clenu N(U) = integral(u*u_x*phi)
    %   F_vec - vektor prave strany

    num_nodes = length(x);

    M     = zeros(num_nodes);
    K     = zeros(num_nodes);
    A_jac = zeros(num_nodes);
    N_vec = zeros(num_nodes, 1);
    F_vec = zeros(num_nodes, 1);

    for e = 1:n
        nodes = (e-1)*p+1 : e*p+1;
        [M_loc, K_loc, A_loc, N_loc, F_loc] = ...
            compute_element_matrices(p, x(nodes), U(nodes), t, params);

        M(nodes, nodes)     = M(nodes, nodes)     + M_loc;
        K(nodes, nodes)     = K(nodes, nodes)     + K_loc;
        A_jac(nodes, nodes) = A_jac(nodes, nodes) + A_loc;
        N_vec(nodes)        = N_vec(nodes)        + N_loc;
        F_vec(nodes)        = F_vec(nodes)        + F_loc;
    end
end

% -------------------------------------------------------------------------

function [M_loc, K_loc, A_loc, N_loc, F_loc] = compute_element_matrices(p, x_elem, U_elem, t, params)

    x_left  = x_elem(1);
    x_right = x_elem(end);
    h       = x_right - x_left;
    Jac     = h/2;

    [gauss_pts, gauss_wts] = gauss_quadrature_info(p+2);

    M_loc = zeros(p+1);
    K_loc = zeros(p+1);
    A_loc = zeros(p+1);
    N_loc = zeros(p+1, 1);
    F_loc = zeros(p+1, 1);

    for i = 1:length(gauss_pts)
        xi = gauss_pts(i);
        w  = gauss_wts(i);

        [phi, dphi_dxi] = lagrange_basis(p, xi);
        dphi_dx = dphi_dxi * (2/h);

        x_phys  = (x_left + x_right)/2 + Jac*xi;
        u_val   = dot(U_elem, phi);
        du_dx   = dot(U_elem, dphi_dx);
        f_val   = params.source(x_phys, t);

        for j = 1:p+1
            N_loc(j) = N_loc(j) + w * u_val * du_dx * phi(j) * Jac;
            F_loc(j) = F_loc(j) + w * f_val * phi(j) * Jac;

            for k = 1:p+1
                M_loc(j,k) = M_loc(j,k) + w * phi(j) * phi(k) * Jac;
                K_loc(j,k) = K_loc(j,k) + w * dphi_dx(j) * dphi_dx(k) * Jac;
                % A_jac = dN/dU: linearizace u*u_x
                % d/dU_k [u*u_x] = phi_k * u_x + u * dphi_dx_k
                A_loc(j,k) = A_loc(j,k) + w * (phi(k)*du_dx + u_val*dphi_dx(k)) * phi(j) * Jac;
            end
        end
    end
end

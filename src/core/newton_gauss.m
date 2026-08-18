function [points, weights] = newton_gauss(num_points)
    % Gaussovy body a vahy pres Newtonovu metodu na Legendreovych polynomech

    tol      = 1e-15;
    max_iter = 50;

    if mod(num_points, 2) == 1
        m        = (num_points - 1) / 2;
        has_zero = true;
    else
        m        = num_points / 2;
        has_zero = false;
    end

    i               = 1:m;
    initial_guesses = cos((i - 0.25) * pi / (num_points + 0.5));
    pos_roots       = zeros(1, m);

    for k = 1:m
        x = initial_guesses(k);

        for iter = 1:max_iter
            [P, dP] = legendre_eval(x, num_points);
            dx      = P / dP;
            x       = x - dx;
            if abs(dx) < tol, break; end
        end

        pos_roots(k) = x;
    end

    pos_roots = sort(pos_roots);

    if has_zero
        points = [-fliplr(pos_roots), 0, pos_roots];
    else
        points = [-fliplr(pos_roots), pos_roots];
    end

    weights = zeros(1, num_points);
    for k = 1:num_points
        [~, dP]    = legendre_eval(points(k), num_points);
        weights(k) = 2 / ((1 - points(k)^2) * dP^2);
    end
end

function [P, dP] = legendre_eval(x, n)
    if n == 0
        P = 1; dP = 0; return;
    elseif n == 1
        P = x; dP = 1; return;
    end

    P_prev  = 1; P_curr  = x;
    dP_prev = 0; dP_curr = 1;

    for j = 1:n-1
        P_next  = ((2*j+1)*x*P_curr  - j*P_prev)  / (j+1);
        dP_next = ((2*j+1)*(P_curr + x*dP_curr) - j*dP_prev) / (j+1);
        P_prev  = P_curr;  P_curr  = P_next;
        dP_prev = dP_curr; dP_curr = dP_next;
    end

    P  = P_curr;
    dP = dP_curr;
end

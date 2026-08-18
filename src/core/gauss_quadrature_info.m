function [points, weights] = gauss_quadrature_info(n)
    % Gauss-Legendre kvadratura na [-1,1]
    % n = pocet bodu - 1  (0=1 bod, 1=2 body, ..., 4=5 bodu)
    % Pro n >= 5 pouziva Newtonovu metodu

    persistent cache_n cache_points cache_weights

    if ~isempty(cache_n) && cache_n == n
        points  = cache_points;
        weights = cache_weights;
        return;
    end

    switch n
        case 0
            points  = 0;
            weights = 2;
        case 1
            points  = [-1/sqrt(3), 1/sqrt(3)];
            weights = [1, 1];
        case 2
            points  = [-sqrt(3/5), 0, sqrt(3/5)];
            weights = [5/9, 8/9, 5/9];
        case 3
            p1 = sqrt((3 - 2*sqrt(6/5))/7);
            p2 = sqrt((3 + 2*sqrt(6/5))/7);
            points  = [-p2, -p1, p1, p2];
            w1 = (18 + sqrt(30))/36;
            w2 = (18 - sqrt(30))/36;
            weights = [w2, w1, w1, w2];
        case 4
            p1 = (1/3)*sqrt(5 - 2*sqrt(10/7));
            p2 = (1/3)*sqrt(5 + 2*sqrt(10/7));
            points  = [-p2, -p1, 0, p1, p2];
            w1 = (322 + 13*sqrt(70))/900;
            w2 = (322 - 13*sqrt(70))/900;
            weights = [w2, w1, 128/225, w1, w2];
        otherwise
            [points, weights] = newton_gauss(n+1);
    end

    cache_n       = n;
    cache_points  = points;
    cache_weights = weights;
end

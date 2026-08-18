function [phi, dphi_dxi] = lagrange_basis(p, xi)

    if p == 1
        phi      = [(1-xi)/2; (1+xi)/2];
        dphi_dxi = [-0.5; 0.5];

    elseif p == 2
        phi      = [0.5*xi*(xi-1); 1-xi^2; 0.5*xi*(xi+1)];
        dphi_dxi = [xi-0.5; -2*xi; xi+0.5];

    elseif p == 3
        phi      = zeros(4,1);
        dphi_dxi = zeros(4,1);

        phi(1) = -9/16 * (xi-1/3) * (xi-1) * (xi+1/3);
        phi(2) =  27/16 * (xi-1)  * (xi-1/3) * (xi+1);
        phi(3) = -27/16 * (xi-1)  * (xi+1/3) * (xi+1);
        phi(4) =   9/16 * (xi-1/3) * (xi+1/3) * (xi+1);

        dphi_dxi(1) = -9/16  * ((xi-1/3)*(xi+1/3) + (xi-1)*(xi+1/3) + (xi-1)*(xi-1/3));
        dphi_dxi(2) =  27/16 * ((xi-1/3)*(xi+1)   + (xi-1)*(xi+1)   + (xi-1)*(xi-1/3));
        dphi_dxi(3) = -27/16 * ((xi+1/3)*(xi+1)   + (xi-1)*(xi+1)   + (xi-1)*(xi+1/3));
        dphi_dxi(4) =   9/16 * ((xi+1/3)*(xi+1)   + (xi-1/3)*(xi+1) + (xi-1/3)*(xi+1/3));

    elseif p == 4
        phi      = zeros(5,1);
        dphi_dxi = zeros(5,1);

        phi(1) =  (2/3) * (xi+0.5) * xi * (xi-0.5) * (xi-1);
        phi(2) = (-8/3) * (xi+1)   * xi * (xi-0.5) * (xi-1);
        phi(3) =  4     * (xi+1)   * (xi+0.5) * (xi-0.5) * (xi-1);
        phi(4) = (-8/3) * (xi+1)   * (xi+0.5) * xi * (xi-1);
        phi(5) =  (2/3) * (xi+1)   * (xi+0.5) * xi * (xi-0.5);

        dphi_dxi(1) = (2/3) * ( (xi)*(xi-0.5)*(xi-1) + (xi+0.5)*(xi-0.5)*(xi-1) + (xi+0.5)*xi*(xi-1) + (xi+0.5)*xi*(xi-0.5) );
        dphi_dxi(2) = (-8/3) * ( (xi)*(xi-0.5)*(xi-1) + (xi+1)*(xi-0.5)*(xi-1) + (xi+1)*xi*(xi-1) + (xi+1)*xi*(xi-0.5) );
        dphi_dxi(3) = 4 * ( (xi+0.5)*(xi-0.5)*(xi-1) + (xi+1)*(xi-0.5)*(xi-1) + (xi+1)*(xi+0.5)*(xi-1) + (xi+1)*(xi+0.5)*(xi-0.5) );
        dphi_dxi(4) = (-8/3) * ( (xi+0.5)*xi*(xi-1) + (xi+1)*xi*(xi-1) + (xi+1)*(xi+0.5)*(xi-1) + (xi+1)*(xi+0.5)*xi );
        dphi_dxi(5) = (2/3) * ( (xi+0.5)*xi*(xi-0.5) + (xi+1)*xi*(xi-0.5) + (xi+1)*(xi+0.5)*(xi-0.5) + (xi+1)*(xi+0.5)*xi );

    else
        error('Stupeň p=%d není implementován (podporováno p=1,2,3,4)', p);
    end
end
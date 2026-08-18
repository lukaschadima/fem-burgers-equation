function T = butcher_tables(method)
% T = butcher_tables(method)
%
% Vraci Butcherovy tabulky pro explicitni, plne implicitni nebo IMEX metody.
%
% Pro explicitni RK (T.type = 'explicit'):
%   T.A  [s x s]  Butcherova matice (strictly lower triangular)
%   T.b  [s x 1]  vahove koeficienty
%   T.c  [s x 1]  casove uzly
%   T.s           pocet stupnu
%
% Pro plne implicitni DIRK (T.type = 'dirk'):
%   T.A  [s x s]  Butcherova matice (lower triangular, A(i,i)>0)
%   T.b  [s x 1]  vahove koeficienty
%   T.c  [s x 1]  casove uzly
%   T.s           pocet stupnu
%   Na kazdem stupni Newton resi cely nelinearni system (konvekce + difuze).
%
% Pro IMEX ARS (T.type = 'imex'), padded tableau:
%   T.Ai, T.bi, T.ci  ... implicitni DIRK cast (difuze, linearni)
%   T.Ae, T.be, T.ce  ... explicitni ERK cast  (konvekce)
%   T.s               ... pocet stupnu (padded)
%
% Dostupna schemata:
%   Explicitni RK   : 'erk1', 'erk2', 'erk3', 'erk4'
%   Plne impl. DIRK : 'dirk2', 'dirk3'
%   IMEX ARS        : 'imex_ars111','imex_ars122','imex_ars232','imex_ars343'
%
% Reference:
%   Ascher, Ruuth, Spiteri (1997), Appl. Numer. Math. [IMEX]
%   Hairer & Wanner (1996), Solving ODEs II           [DIRK]

switch lower(method)

    % =========================================================
    %  EXPLICITNI RUNGE-KUTTA
    % =========================================================

    case 'erk1'   % Forward Euler, rad 1
        T.A = 0;
        T.b = 1;
        T.c = 0;
        T.s = 1;
        T.type = 'explicit';

    case 'erk2'   % Heunova metoda (explicitni trapezoid), rad 2
        T.A = [0 0;
               1 0];
        T.b = [1/2; 1/2];
        T.c = [0; 1];
        T.s = 2;
        T.type = 'explicit';

    case 'erk3'   % Klasicke RK3 (Kuttova metoda), rad 3
        T.A = [0    0   0;
               1/2  0   0;
              -1    2   0];
        T.b = [1/6; 2/3; 1/6];
        T.c = [0; 1/2; 1];
        T.s = 3;
        T.type = 'explicit';

    case 'erk4'   % Klasicke RK4, rad 4
        T.A = [0    0    0   0;
               1/2  0    0   0;
               0    1/2  0   0;
               0    0    1   0];
        T.b = [1/6; 1/3; 1/3; 1/6];
        T.c = [0; 1/2; 1/2; 1];
        T.s = 4;
        T.type = 'explicit';

    % =========================================================
    %  PLNE IMPLICITNI SDIRK SCHEMATA
    %  Na kazdem stupni i resime Newtonovou metodou:
    %    R(U_i) = M*U_i + dt*A(i,i)*(N(U_i)+nu*K*U_i) - rhs_i = 0
    %  kde rhs_i obsahuje explicitni prispevky z predchozich stupnu.
    %  Bezpodminecna stabilita, ale Newton = drahy krok.
    % =========================================================

    case 'dirk2'
        % SDIRK2: 2-stage, 2nd order, L-stabilni
        % Alexander (1977), Hairer & Wanner (1996) p.100
        % gamma = (2 - sqrt(2)) / 2  ~  0.2929
        gam = (2 - sqrt(2)) / 2;
        T.A = [gam,       0  ;
               1 - gam,   gam];
        T.b = [1 - gam; gam];
        T.c = [gam; 1];
        T.s = 2;
        T.type = 'dirk';

    case 'dirk3'
        % SDIRK3: 3-stage, 3rd order, L-stabilni
        % Stejny gamma jako v ARS(3,4,3) -- stridne-stabilni DIRK cast.
        % Hairer & Wanner (1996), Tab. IV.6.5, s.100
        % gamma = stredni koren 6x^3 - 18x^2 + 9x - 1 = 0
        gam = 0.4358665215484241;
        b1  = -1.5*gam^2 + 4*gam - 0.25;    %  ~1.2085
        b2  =  1.5*gam^2 - 5*gam + 1.25;    %  ~-0.6444

        T.A = [gam,          0,     0;
               (1-gam)/2,    gam,   0;
               b1,           b2,    gam];
        T.b = [b1; b2; gam];   % stiffly accurate: b = A(end,:)
        T.c = [gam; (1+gam)/2; 1];
        T.s = 3;
        T.type = 'dirk';

    % =========================================================
    %  IMEX-ARS SCHEMATA  (Ascher, Ruuth, Spiteri 1997)
    %  Padded tableau: prvni stupen trivialní (Ai(1,:)=0, Ae(1,:)=0).
    %  Implicitni cast: difuze -> linearni system (M + dt*aii*nu*K).
    %  Explicitni cast: konvekce -> bez Newtona.
    % =========================================================

    case 'imex_ars111'
        % ARS (1,1,1): Forward-Backward Euler, rad 1
        % Ekvivalentni se stavajicim step_imex_euler.
        T.Ai = [0  0;
                0  1];
        T.bi = [0; 1];
        T.ci = [0; 1];
        T.Ae = [0  0;
                1  0];
        T.be = [1; 0];
        T.ce = [0; 1];
        T.s    = 2;
        T.type = 'imex';

    case 'imex_ars122'
        % ARS (1,2,2): IMEX Midpoint, rad 2
        T.Ai = [0    0  ;
                0    1/2];
        T.bi = [0; 1];
        T.ci = [0; 1/2];
        T.Ae = [0    0;
                1/2  0];
        T.be = [0; 1];
        T.ce = [0; 1/2];
        T.s    = 2;
        T.type = 'imex';

    case 'imex_ars232'
        % ARS (2,3,2): L-stabilni DIRK, rad 2
        gam = (2 - sqrt(2)) / 2;
        del = -2*sqrt(2) / 3;
        T.Ai = [0,     0,     0;
                0,     gam,   0;
                0,   1-gam,   gam];
        T.bi = [0; 1-gam; gam];
        T.ci = [0; gam; 1];
        T.Ae = [0,    0,      0;
                gam,  0,      0;
                del,  1-del,  0];
        T.be = [0; 1-gam; gam];
        T.ce = [0; gam; 1];
        T.s    = 3;
        T.type = 'imex';

    case 'imex_ars343'
        % ARS (3,4,3): L-stabilni, rad 3
        gam = 0.4358665215484241;
        b1  = -1.5*gam^2 + 4*gam - 0.25;
        b2  =  1.5*gam^2 - 5*gam + 1.25;
        T.Ai = [0,           0,          0,    0;
                0,           gam,        0,    0;
                0,  (1-gam)/2,            gam,  0;
                0,           b1,         b2,   gam];
        T.bi = [0; b1; b2; gam];
        T.ci = [0; gam; (1+gam)/2; 1];
        a31 =  0.3212788860;
        a32 =  0.3966543747;
        a41 = -0.105858296;
        a42 =  0.5529291479;
        a43 =  0.5529291479;
        T.Ae = [0,    0,    0,    0;
                gam,  0,    0,    0;
                a31,  a32,  0,    0;
                a41,  a42,  a43,  0];
        T.be = [0; b1; b2; gam];
        T.ce = [0; gam; (1+gam)/2; 1];
        T.s    = 4;
        T.type = 'imex';

    otherwise
        error(['butcher_tables: nezname schema "%s".\n' ...
               'Dostupna: erk1, erk2, erk3, erk4, dirk2, dirk3,\n' ...
               '          imex_ars111, imex_ars122, imex_ars232, imex_ars343'], method);
end
end

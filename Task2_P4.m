q10 = 0;
q20 = pi;

Ts = 0.001;

X0=[0; 0.95*pi; 0; 0];

M = [Ph(1)  Ph(6)*cos(q20);
            Ph(6)*cos(q20)  Ph(2)];

F = [Ph(3)  0;
              0          Ph(4)];

K =  [0  0;
      0  Ph(5)*cos(q20)];

G = [1; 0];


A = [zeros(2,2)  eye(2);
     -M\K   -M\F];

B= [zeros(2,1);
    M\G];
 
C=[1 0 0 0;
   0 1 0 0];

p = [-10+10*j -10-10*j -15 -18];

alpha = poly(p);
p_op= eig(A);
a=poly(p_op);

Kc  = [ alpha(5) - a(5), ...
        alpha(4) - a(4), ...
        alpha(3) - a(3), ...
        alpha(2) - a(2) ];

Pc  = [ B, A*B, A^2*B, A^3*B ];

Bc  = [0; 0; 0; 1];

Ac  = [ 0    1    0    0;
        0    0    1    0;
        0    0    0    1;
       -a(5) -a(4) -a(3) -a(2) ];

Pcb = [ Bc, Ac*Bc, Ac^2*Bc, Ac^3*Bc ];

T_1 = Pcb / Pc;      % = Pcb*inv(Pc) = T^{-1}

Kf  = Kc * T_1;



Lo=place(A',C', 6*p);
Lo=Lo';

Po    = obsv(A, C);      % observability matrix
r_Po  = rank(Po);        % rank of Po
n     = size(A,1);


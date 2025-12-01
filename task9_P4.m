A11 = zeros(2);          % 0_{2x2}
A12 = eye(2);            % I_{2x2}
A21 = -M\K;
A22 = -M\F;

B1 = zeros(2,1);
B2 = M\G;


p_red = [-40 -50];


Lro = place(A22', A12', p_red).'; 

Aro = A22 - Lro*A12;                                    
Bro = [ B2 - Lro*B1 , ...                               
        A22*Lro - Lro*A12*Lro + A21 - Lro*A11 ];        

Cro = eye(2);                                           
Dro = [ zeros(2,1) , Lro ];                             



z0 = -Lro * X0;

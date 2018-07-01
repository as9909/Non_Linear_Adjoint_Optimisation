syms x y z dt n1 n2 l1 l2 l3
% a3*dt=n1, (a2+a3)*dt=n2,  l1=u1-u2, l2= u1-u3, l3=u1-u4, a=x, b=y, c=z
eqns = [x*n1^3+y*n1^2+z*n1+l1==0, ...
        x*n2^3+y*n2^2+z*n2+l2==0, ...
        x*dt^3+y*dt^2+z*dt+l3==0];

S = solve(eqns);
sol = collect(simplify(expand([S.x; S.y; S.z])));

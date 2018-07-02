syms x y z dt n1 n2 l1 l2 l3
% a3*dt=n1, (a2+a3)*dt=n2,  l1=u1-u2, l2= u1-u3, l3=u1-u4, a=x, b=y, c=z
eqns = [x*n1^3+y*n1^2+z*n1+l1==0, ...
        x*n2^3+y*n2^2+z*n2+l2==0, ...
        x*dt^3+y*dt^2+z*dt+l3==0];

S = solve(eqns);
x_sol=(simplify((S.x)));
y_sol=(simplify((S.y)));
z_sol=(simplify((S.z)));


% The answers have been observed and put into condensed form manually
beta=dt*n1*n2*(dt-n1)*(dt-n2)*(n1-n2);

x1_sol=( l1*n2*dt*(dt-n2)-l2*n1*dt*(dt-n1)+l3*n1*n2*(n2-n1) )/beta;

y1_sol=( l2*n1*dt*(dt^2-n1^2)-l1*n2*dt*(dt^2-n2^2)+l3*n1*n2*(n1^2-n2^2) )/beta;

z1_sol=( l1*n2^2*dt^2*(dt-n2)-l2*n1^2*dt^2*(dt-n1)+l3*n1^2*n2^2*(n2-n1) )/beta;

simplify(expand(x_sol-x1_sol))

simplify(expand(y_sol-y1_sol))

simplify(expand(z_sol-z1_sol))
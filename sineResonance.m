function out = sineResonance(frequency) %first using eta = 0.1, let A, F^2 = 1

syms x(t)
syms t
syms f(t)
f(t) = sin(frequency*t);


eqn  = diff(x,t,2) + 0.1*diff(x,t) + x == f(t);

%we want to fix energy
%it has been calculated that in order to fix energy (ignoring a 
%transcendental contribution, should use T = 2*F^2/A^2 where F is the total
%energy of the beam. For now, we use F^2 = A^2 = 1, and are interested in
%the resonance curve of this system



%EDIT: DO NOT IGNORE TRANSCENDENTAL COMPONENT
syms Tsqrt;
allowedPds = vpasolve((Tsqrt^2/2 - sin(2*pi*Tsqrt^2*frequency)/(4*frequency)) == 1,Tsqrt);

actualpd = double(allowedPds^2)



[V] = odeToVectorField(eqn);
M = matlabFunction(V, 'vars', {'t', 'Y'});
sol = ode45(M, [0 actualpd], [0 0]);

fplot(@(x)deval(sol, x, 1), [0 actualpd], '-b')

sz = size(sol.y, 2);


finalPos = sol.y(1, sz);
finalVel = sol.y(2, sz);



out = [finalPos.^2 + finalVel.^2]
end

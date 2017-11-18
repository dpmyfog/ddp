function out = sineArray(eta, omega, A) %a is the amplitude of the forcing function  
%takes in eta, omega and A and returns the amplitude of the next
%oscillation

%start with some solution
syms t;
syms x(t);


syms Tsqrt;
allowedPds = vpasolve((Tsqrt^2/2 - sin(2*pi*Tsqrt^2*frequency)/(4*frequency)) == 1,Tsqrt)

actualpd = double(allowedPds^2)

x(t) = (A*exp(-t*(eta/2 - ((eta - 2)*(eta + 2))^(1/2)/2))*exp((eta*t)/2 - (t*(eta^2 - 4)^(1/2))/2)*(sin(omega*t)*(eta/2 - (eta^2 - 4)^(1/2)/2) - omega*cos(omega*t)))/((eta^2 - 4)^(1/2)*((eta/2 - (eta^2 - 4)^(1/2)/2)^2 + omega^2)) - (4*A*omega*exp(-t*(eta/2 + ((eta - 2)*(eta + 2))^(1/2)/2)))/((eta^2 - 4)^(1/2)*(2*eta*(eta^2 - 4)^(1/2) + 2*eta^2 + 4*omega^2 - 4)) - (4*A*omega*exp(-t*(eta/2 - ((eta - 2)*(eta + 2))^(1/2)/2)))/((eta^2 - 4)^(1/2)*(2*eta*(eta^2 - 4)^(1/2) - 2*eta^2 - 4*omega^2 + 4)) - (A*exp(-t*(eta/2 + ((eta - 2)*(eta + 2))^(1/2)/2))*exp((eta*t)/2 + (t*(eta^2 - 4)^(1/2))/2)*(sin(omega*t)*(eta/2 + (eta^2 - 4)^(1/2)/2) - omega*cos(omega*t)))/((eta^2 - 4)^(1/2)*((eta/2 + (eta^2 - 4)^(1/2)/2)^2 + omega^2));
xprime(t) = diff(x, t);

%now need to integrate 

powerLost = double(int(xprime(t)*A*eta*sin(omega*t), 0, actualpd));
coeff = -4*omega*eta/(sin(4*omega) - 4*omega);

out = sqrt(A^2 - coeff*powerLost);

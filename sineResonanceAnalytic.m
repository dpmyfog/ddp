function out = sineResonanceAnalytic(frequency) %first using eta = 0.1, let A, F^2 = 1

syms t
syms x(t)
eta = 0.1;


x(t)  = (exp(-t*(eta/2 - ((eta - 2)*(eta + 2))^(1/2)/2))*exp((eta*t)/2 - (t*(eta^2 - 4)^(1/2))/2)*(sin(frequency*t)*(eta/2 - (eta^2 - 4)^(1/2)/2) - frequency*cos(frequency*t)))/((eta^2 - 4)^(1/2)*((eta/2 - (eta^2 - 4)^(1/2)/2)^2 + frequency^2)) - (4*frequency*exp(-t*(eta/2 + ((eta - 2)*(eta + 2))^(1/2)/2)))/((eta^2 - 4)^(1/2)*(2*eta*(eta^2 - 4)^(1/2) + 2*eta^2 + 4*frequency^2 - 4)) - (4*frequency*exp(-t*(eta/2 - ((eta - 2)*(eta + 2))^(1/2)/2)))/((eta^2 - 4)^(1/2)*(2*eta*(eta^2 - 4)^(1/2) - 2*eta^2 - 4*frequency^2 + 4)) - (exp(-t*(eta/2 + ((eta - 2)*(eta + 2))^(1/2)/2))*exp((eta*t)/2 + (t*(eta^2 - 4)^(1/2))/2)*(sin(frequency*t)*(eta/2 + (eta^2 - 4)^(1/2)/2) - frequency*cos(frequency*t)))/((eta^2 - 4)^(1/2)*((eta/2 + (eta^2 - 4)^(1/2)/2)^2 + frequency^2));

%we want to fix energy
%it has been calculated that in order to fix energy (ignoring a
%transcendental contribution, should use T = 2*F^2/A^2 where F is the total
%energy of the beam. For now, we use F^2 = A^2 = 1, and are interested in
%the resonance curve of this system

syms Tsqrt;
allowedPds = vpasolve((Tsqrt^2/2 - sin(2*pi*Tsqrt^2*frequency)/(4*frequency)) == 1,Tsqrt);

actualpd = double(allowedPds^2)

xprime(t) = diff(x, t)

x(t) = rewrite(x, 'sincos')
xprime(t) = rewrite(xprime, 'sincos')

f = @(t) x(t)
fprime = @(t) xprime(t)
finalPos = f(actualpd)
finalVel = fprime(actualpd)

RePos = simplify(collect(real(finalPos)))
ReVel = simplify(collect(real(finalVel)))

imPos = simplify(collect(imag(finalPos)))
imVel = simplify(collect(imag(finalVel)))

out = double(RePos).^2 + double(ReVel).^2;
end

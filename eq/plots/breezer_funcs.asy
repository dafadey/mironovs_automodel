real ch(real x) {return (exp(x) + exp(-x)) * .5;}

real sh(real x) {return (exp(x) - exp(-x)) * .5;}

real breezer(real x, real z, real a, real d)
{
	real E = x - (a^2 - 3 * d^2) * z;
	real P = x - (3 * a^2 - d^2) * z;
	
	return 2*a*d*(a*sh(a*E)*sin(d*P)-d*ch(a*E)*cos(d*P))/(a^2*sin(d*P)^2+d^2*ch(a*E)^2);
}

real breezerPhase(real x, real phase, real a, real d)
{
	real E = x;
	real P = x - phase;

	// E = x - (a^2-3*d^2)*z
	// P = x - (3a^2-d^2)*z
	
	return 2*a*d*(a*sh(a*E)*sin(d*P)-d*ch(a*E)*cos(d*P))/(a^2*sin(d*P)^2+d^2*ch(a*E)^2);
}



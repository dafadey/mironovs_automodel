import breezer_funcs;
import graph;

picture p;

guide ge;

real d=230;
real a=100;
real phase=0;
for(real par=70; par < 150; par+=3)
//for(real par=0; par < 2*pi/d; par+=0.1/d)
//for(real par=150; par < 300; par+=10)
{
  guide g;
  real W=0;
  real dx = 0.1/d;
  for(real x=-6*pi/d; x < 6*pi/d; x+=dx)
  {
		real v = breezerPhase(x, 0, par, d);
		W += v*v*dx;
		g = g -- (x, v);
	}
	ge = ge -- (par, W);
	draw(p,g);
}

xaxis(p, "$\eta$", BottomTop, LeftTicks);
yaxis(p, "$u$", LeftRight, RightTicks);
p = shift(-point(p, SW)) * p;
size(p, 7cm, 6cm, point(p,SW), point(p,NE));
add(p.fit());

picture pe;
draw(pe, ge);
xaxis(pe, "parameter", BottomTop, LeftTicks);
yaxis(pe, "$W$", LeftRight, RightTicks);
pe = shift(-point(pe, SW)) * pe;
size(pe, 7cm, 6cm, point(pe,SW), point(pe,NE));
add(pe.fit(), (8.7cm,0));

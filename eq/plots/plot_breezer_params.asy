import graph;
import breezer_funcs;

real q=4e-7;

real getW(real a, real d, real phase, real z0)
{
  real W=0;
  real dx = 0.0001;
  for(real x=-0.2; x<0.2; x+=dx)
    W += (q^0.5 * z0 * breezerPhase(x, phase, a, d))^2 * dx;
  return W;
}

file inf = input("breezer_params.dat");

real[] raw = inf;

int Nfields=6;

int Nrec = floor(raw.length/Nfields);

real[] zs;
real[] as;
real[] ds;
real[] ps;
real[] x0s;
real[] fs;


for(int i=0; i<Nrec; ++i)
{
	zs.push(raw[i*Nfields]);
	as.push(raw[i*Nfields+1]);
	ds.push(raw[i*Nfields+2]);
	ps.push(raw[i*Nfields+3]);
	x0s.push(raw[i*Nfields+4]);
	fs.push(raw[i*Nfields+5]);
}

guide ags;
guide dgs;
guide pgs;
guide x0gs;
guide fgs;

guide avcentergs;
guide avphasegs;

guide awgs;

guide vcentergs;
guide vphasegs;

real pp=0;
real p0=0;
real dp=0;

for(int i=0; i<Nrec; ++i) {
	real z = zs[i];
	ags=ags--(zs[i],as[i]);
	dgs=dgs--(zs[i],ds[i]);
	real d = ds[i];
  real p = ps[i];//*d;
	//p = p - floor(p/pi)*pi;
	//if(p < pp)
	//	p0+=pi;
  if(i>0)
  {
    real vc = (x0s[i]-x0s[i-1])/(zs[i]-zs[i-1]);
    real vph = (p-pp)/(zs[i]-zs[i-1]);//(ps[i+1]-ps[i])/(zs[i+1]-zs[i]);
    vcentergs = vcentergs -- ((zs[i]+zs[i-1])*0.5, vc);
    if(abs(ps[i]-ps[i-1])*ds[i] < pi/10.)
    vphasegs = vphasegs -- ((zs[i]+zs[i-1])*0.5, vc + vph);
  }
	pp=p;
  dp=d;
	pgs=pgs--(zs[i],(p+p0)/ds[i]);
	
	x0gs=x0gs--(zs[i],x0s[i]);
	fgs=fgs--(zs[i],abs(fs[i]));
  avphasegs = avphasegs -- (zs[i], (3*as[i]^2-ds[i]^2)*q);
  avcentergs = avcentergs -- (zs[i], (as[i]^2-3*ds[i]^2)*q);

  awgs = awgs -- (zs[i], getW(as[i], ds[i], ps[i], zs[i]));
}

void add_pic(guide[] garr, pair pos, string name, pen[] pens = {currentpen}) {
  picture p;
  for(int i=0; i<garr.length; ++i)
    draw(p,garr[i], pens[i%pens.length],marker(unitcircle, pens[i%pens.length]));
  xaxis(p, "$z$", BottomTop, LeftTicks);
	yaxis(p, name, LeftRight, RightTicks);
	p = shift(-point(p, SW)) * p;
  size(p, 5cm, 4cm, point(p, SW), point(p, NE));
  add(p.fit(), pos);
}

add_pic(new guide[]{ags, dgs}, (0cm,0cm), "$a$, $\delta$", new pen[] {black, red});
add_pic(new guide[] {pgs}, (14cm,0), "$\phi$");
add_pic(new guide[] {x0gs}, (0cm,6cm), "$x_0$");

add_pic(new guide[] {vcentergs, avcentergs, vphasegs, avphasegs}, (14cm,6cm), "$\displaystyle \left. \frac{\partial\eta}{\partial z} \right\vert_E$, $\displaystyle \left. \frac{\partial\eta}{\partial z} \right\vert_P$", new pen[] {black, red, blue, 0.7*green});

add_pic(new guide[]{awgs}, (0cm,-6cm), "W");



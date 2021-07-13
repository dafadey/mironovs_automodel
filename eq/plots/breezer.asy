import graph;
import breezer_funcs;
default(fontsize(8pt));
file fin=input("signal.dat");
real[] raw = fin;
write("values count is ", raw.length);
pair[][] signals;
int Nx;
int i=0;
real[] zs;
int Nsigs=floor(raw[i]);
i+=1;
write("N signals ", Nsigs);
{
  for(int j=0;j<Nsigs;++j)
  {
    Nx=floor(raw[i]);
    i+=1;
    zs.push(raw[i]);
    i+=1;
    pair[] signal;
    for(int ii=0; ii<Nx; ++ii)
    {
      signal.push((raw[i], raw[i+1]));
      i+=2;	
    }
    signals.push(signal);
  }
}

int limit = 42;

real getControlPointAnalytic(real _xmax, real _dx, real phase, real a, real d)
{
	real xmax=_xmax;
	real dx = _dx;
	while(dx > pi/a/100)
	{
		real x0=xmax-dx;
		real x1=xmax+dx;
		dx=(x1 - x0)*0.125;
		for(real x = x0; x<x1; x+=dx)
		{
			if(abs(breezerPhase(x, phase, a, d))>abs(breezerPhase(x+dx, phase, a, d)) && 
		   	   abs(breezerPhase(x, phase, a, d))>abs(breezerPhase(x-dx, phase, a, d)))
			{
				xmax = x;
				break;
			}
		}
	}
	return xmax;
}

pair[] getControlPointsAnalytic(real phase, real a, real d)
{
  //write("getControlPointsAnalytic a="+string(a)+" d="+string(d)+" phi="+string(phase));
	real[] xmaxes;
	real dx = pi/d*0.01;//0.128;
	for(real x=-pi/d*3; x<pi/d*3; x+=dx)
	{
		if(abs(breezerPhase(x, phase, a, d))>abs(breezerPhase(x+dx, phase, a, d)) && 
	   	   abs(breezerPhase(x, phase, a, d))>abs(breezerPhase(x-dx, phase, a, d)))
			xmaxes.push(x);
	}

	real absmax=0;
	int imax=0;
	for(int i=0; i<xmaxes.length; ++i)
	{
		real va = abs(breezerPhase(xmaxes[i], phase, a, d));
		if(absmax < va)
		{
			absmax = va;
			imax = i;
		}
	}

	pair[] res;
	if(imax>0 && imax+1<xmaxes.length)
	{
		for(int i=imax-1; i<=imax+1; ++i) {
			real xmax=getControlPointAnalytic(xmaxes[i], dx, phase, a, d);
			res.push((xmax, breezerPhase(xmax, phase, a, d)));
		}
	}	
	return res;
}

int getMaxPos(pair[] signal)
{
  real absmax=0;
  int maxpos=0;
  for(int i=0;i<Nx;++i)
  {
    if(absmax < abs(signal[i].y))
    {
      absmax=abs(signal[i].y);
      maxpos = i;
    }
  }
  return maxpos;
}

real[] getPeriod(pair[] signal)
{
  real xmin=signal[0].x;
  real xmax=signal[signal.length-1].x;
  
  write("xmin=",xmin);
  write("xmax=",xmax);
  write("Nx=", Nx);
  
  int maxpos = getMaxPos(signal);
  real[] res;
  int i1,i0;
  for(int i=maxpos; i>0; --i)
  {
    if(signal[i-1].y*signal[i].y<=0)
    {
      i0=i;
      break;
    }
  }
  real x0 = (i0 - abs(signal[i0].y)/(abs(signal[i0].y)+abs(signal[i0-1].y)))/Nx*(xmax-xmin);
  for(int i=maxpos; i<Nx-1; ++i)
  {
    if(signal[i+1].y*signal[i].y<=0)
    {
      i1=i;
      break;
    }
  }
  real x1 = (i1 + abs(signal[i1].y)/(abs(signal[i1].y)+abs(signal[i1+1].y)))/Nx*(xmax-xmin);
  res.push(x0);
  res.push(x1);
  res.push(maxpos/Nx*(xmax-xmin));
  return res;
}

pair max3(pair[] signal, int i) {
  //y=a*x^2+b*x+c
  //y0 = a*dx*2-b*dx+c
  //y1 = c
  //y2 = a*dx*2+b*dx+c
  //y2-y0 = 2*b*dx
  //y2+y0-2*y1 = 2*a*dx^2
  
  real dx = signal[1].x-signal[0].x;
  real c = signal[i].y;
  real b = (signal[i+1].y-signal[i-1].y)/2/dx;
  real a = (signal[i+1].y+signal[i-1].y-2*signal[i].y)/2/dx^2;
  
  //2*a*xmax+b=0;
  real xmax = signal[i].x-b/2/a;
  real ymax = -b^2/4/a+c;
  return (xmax, ymax);
}

pair[] getControlPoints(pair[] signal)
{
  real dx = signal[1].x-signal[0].x;
  int[] xmaxes;
  real[] zeros = getPeriod(signal);
  write("period=", zeros[1]-zeros[0]);
  real d = pi/(zeros[1]-zeros[0]);
  int dxi = floor(pi/d*3 / dx);
  int imax = floor(zeros[2]/dx);
  write("dx=",dx);	
  write("dxi=",dxi);	
  write("imax=",imax);	
  for(int i=imax-dxi; i<imax+dxi; i+=1)
  {
    if(abs(signal[i])>abs(signal[i+1]) && 
         abs(signal[i])>abs(signal[i-1]))
      xmaxes.push(i);
  }
  write("&&&");
  for(int imax : xmaxes)
		write(string(imax) + " : " + string(signal[imax].x)+" "+string(signal[imax].y));

  real absmax=0;
  int maxid=0;
  for(int i=0; i<xmaxes.length; ++i)
  {
    real va = abs(signal[xmaxes[i]]);
    if(absmax < va)
    {
      absmax = va;
      maxid = i;
    }
  }

	write("maxid=",maxid);
  pair[] res;

	for(int i=max(maxid-1,0); i<=min(maxid+1,xmaxes.length-1); ++i)
		res.push(max3(signal, xmaxes[i]));

  return res;
}

real q = 4e-7;
real beta;
real a,d,phase;

real getYAbsMax(pair[] arr)
{
	real res=0;
	for(pair v : arr)
		res=max(abs(v.y), res);
	return res;
}

real getXPeriod3(pair[] arr)
{
  //write("getXPeriod3 arr.length=", arr.length);
	if(arr.length >= 3)
		return (abs(arr[1].x-arr[0].x)+abs(arr[2].x-arr[1].x)) * 0.5;
  else if(arr.length == 2)
		return abs(arr[1].x-arr[0].x);
	else
		return 0;
}

real getXPeriod2(pair[] arr)
{
  /*
	if(arr.length == 3)
    write("getXPeriod2 arr.length="+string(arr.length)+" choosing "+ (abs(arr[0].y) > abs(arr[2].y) ? "1 and 0" : "2 and 1"));
	if(arr.length == 2)
    write("getXPeriod2 arr.length="+string(arr.length)+" choosing 1 and 0");
	*/
  if(arr.length == 3)
		return abs(arr[0].y) > abs(arr[2].y) ?  abs(arr[1].x-arr[0].x) : abs(arr[2].x-arr[1].x);
  else if(arr.length == 2)
		return abs(arr[1].x-arr[0].x);
	else
		return 0;
}

real getMaxRatio2(pair[] arr)
{
  if(arr.length==3)
    return abs(arr[0].y) > abs(arr[2].y) ? abs(arr[1].y)/abs(arr[0].y) : abs(arr[2].y)/abs(arr[1].y);
  else if(arr.length==2)
    return abs(arr[1].y)/abs(arr[0].y);
  else
    return 1;
}

real getMaxRatio3(pair[] arr)
{
  return abs(arr[2].y)/abs(arr[0].y);
}

bool findBreezerFor(pair[] signal, real zs, bool fake_run=false)
{
	beta = q^0.5 * zs;

  int maxpos=getMaxPos(signal);

	write("maxpos=",maxpos);

  //get period

  pair[] maxes=getControlPoints(signal);
  write("***");
  write(maxes);
	
	for(int i=0;i<maxes.length; ++i)
		maxes[i] = (maxes[i].x, maxes[i].y / beta);	

	if(maxes.length < 2)
		return false;
	
  d = pi/getXPeriod3(maxes);
		
  write("***d=",d);
  a=getYAbsMax(maxes)/4;
 
  write("***a=",a);
  phase=0;

  pair[] maxesA=getControlPointsAnalytic(phase,a,d);
	write("***");

  pair[] normalize(pair[] res)
  {
    pair[] res_target;
    for(int i=0;i<res.length;++i)
      res_target.push((res[i].x,res[i].y/res[1].y));
    return res_target;
  }

  pair[] target = normalize(maxes);
	if(target.length<2)
		return false;

  write("-------->");
  
  if(fake_run)
    return true;

  int nmaxes;
  real maxRatio;
  real maxXPeriod;
  if(maxes.length == 3)
  {
    nmaxes = 3;
    maxRatio = getMaxRatio3(maxes);
    maxXPeriod = getXPeriod3(maxes);
  }
  else if(maxes.length == 2)
  {
    nmaxes = 2;
    maxRatio = getMaxRatio2(maxes);
    maxXPeriod = getXPeriod2(maxes);
  }
  else
    return false;

  
  bool solve_a()
  {
		int iters=0;
    real da = a/4;
    real dir_a=1;
    while(true)
    {
      //write("a=",a);
			if(iters>limit)
				return false;
			pair[] maxesA=getControlPointsAnalytic(phase,a,d);
			if(maxesA.length<2)
				return false;
      a=a+dir_a*da;
      pair[] maxesA_a=getControlPointsAnalytic(phase,a,d);
			if(maxesA_a.length<2) {
				dir_a=-dir_a;
				iters+=1;
				continue;
			}
      //write("A=", getYAbsMax(maxesA));
      //if(abs(min(abs(maxesA[0].y),abs(maxesA[2].y))-min(abs(target[0].y),abs(target[2].y))) < abs(min(abs(maxesA_a[0].y),abs(maxesA_a[2].y))-min(abs(target[0].y),abs(target[2].y))))
      if(abs(getYAbsMax(maxesA)- getYAbsMax(maxes)) < abs(getYAbsMax(maxesA_a)- getYAbsMax(maxes)))
      {
        dir_a=-dir_a;
        da=da/2;
      }
      if(da < 0.3)
        break;

			iters+=1;
    }
    
    return true;
  }

  bool solve_d()
  {
		//write("d=",d);
		int iters=0;
    real dd = d/4;
    real dir_d=1;
    //write("solving d a="+string(a)+" d="+string(d)+" phase="+string(phase));
    while(true)
    {
			if(iters>limit)
				return false;
      pair[] maxesA=normalize(getControlPointsAnalytic(phase,a,d));
			if(maxesA.length<2)
				return false;
      d=d+dir_d*dd;
      pair[] maxesA_d=normalize(getControlPointsAnalytic(phase,a,d));
			if(maxesA_d.length<2) {
				dir_d=-dir_d;
				iters+=1;
				continue;
			}
			
      real criteriaA = (nmaxes == 2 || maxesA.length==2) ? abs(getXPeriod2(maxesA)-maxXPeriod) : abs(getXPeriod3(maxesA)-maxXPeriod);
      real criteriaA_d = (nmaxes == 2 || maxesA_d.length==2) ? abs(getXPeriod2(maxesA_d)-maxXPeriod) : abs(getXPeriod3(maxesA_d)-maxXPeriod);
      /*
      write("----------->");
      write("-->d=",d);
      write("-->crA=",criteriaA);
      write("-->crA_d=",criteriaA_d);
      write("-->dd=",dd);
      write("--->getXPeriod2(maxesA_d)="+string(getXPeriod2(maxesA_d))+" "+string(maxXPeriod));
      write("--->getXPeriod2(maxesA)="+string(getXPeriod2(maxesA))+" "+string(maxXPeriod));
      write("");
      */
			if(criteriaA < criteriaA_d || abs(d+dir_d*dd)<0.1*abs(d))
      {
        dir_d=-dir_d;
        dd=dd/2;
      }

      if(dd < 0.3)
        break;

			iters+=1;
    }
    return true;
  }

  bool solve_p()
  {
		//write("phase=",phase);
		int iters=0;
    real dp = 0.25/d;
    real dir_p=1;
        
    while(true)
    {
			if(iters>limit)
				return false;
      pair[] maxesA=normalize(getControlPointsAnalytic(phase,a,d));
      if(maxesA.length<2)
				return false;
      phase=phase+dir_p*dp;
      pair[] maxesA_p=normalize(getControlPointsAnalytic(phase,a,d));
      if(maxesA_p.length<2) {
				dir_p=-dir_p;
				phase=phase+dir_p*dp;
				iters+=1;
				continue;
			}
			
			//real criteriaA = (maxesA.length==2 || target.length==2) ? abs(maxesA[1].y/maxesA[0].y-target[1].y/target[0].y) : abs(maxesA[2].y/maxesA[0].y-target[2].y/target[0].y);
			
			//real criteriaA_p = (maxesA_p.length==2 || target.length==2) ? abs(maxesA_p[1].y/maxesA_p[0].y-target[1].y/target[0].y) : abs(maxesA_p[2].y/maxesA_p[0].y-target[2].y/target[0].y);
			
      real criteriaA = (nmaxes == 2 || maxesA.length==2) ? abs(getMaxRatio2(maxesA)-maxRatio) : abs(getMaxRatio3(maxesA)-maxRatio);
      real criteriaA_p = (nmaxes == 2 || maxesA_p.length==2) ? abs(getMaxRatio2(maxesA_p)-maxRatio) : abs(getMaxRatio3(maxesA_p)-maxRatio);
      
      if(criteriaA < criteriaA_p)
      {
        dir_p=-dir_p;
        dp=dp/2;
      }
      if(dp < 0.001/d)
        break;

			iters+=1;
    }
    return true;
  }

  real critd;
  real critp;
  real crita;

  bool solve()
  {
		real a_p = a;
		real d_p = d;
		int iters=0;
    while(true)
    {
			//write(iters);
			if(iters>limit)
				return true;//false;

      //d_p = d;
      if(!solve_d())
				a = (a+a_p)*0.5;

			solve_p();
      //if(!solve_p())
      //  return true;//false;

			a_p = a;
			solve_a();
      //if(!solve_a())
			//	return true;//false;

      pair[] maxesA=getControlPointsAnalytic(phase,a,d);
      
      //write("*A=", getYAbsMax(maxesA));
      
      if(maxesA.length<2) {
				write("bad a:", a);
				write("bad d:", d);
				write("bad phase:", phase);
        d = (d+d_p)*0.5;
				write("maxesA.length<2");
				iters+=1;
        //continue;
				return false;
			}
      
      critd =(nmaxes == 2 || maxesA.length==2) ? abs(getXPeriod2(maxesA)-maxXPeriod) : abs(getXPeriod3(maxesA)-maxXPeriod);
      critp = (nmaxes == 2 || maxesA.length==2) ? abs(getMaxRatio2(maxesA)-maxRatio) : abs(getMaxRatio3(maxesA)-maxRatio);
      crita =abs(getYAbsMax(maxesA)-getYAbsMax(maxes));
      write("iter# "+string(iters)+" cra="+string(crita)+" crd="+string(critd)+" crp="+string(critp));
      
      //write("aright: ", abs(maxesA[2].y-target[2].y));
      //write("aleft: ", abs(maxesA[0].y-target[0].y));
      //write("phase: ", abs(maxesA[2].x-maxesA[0].x-(target[2].x-target[0].x)));
      if(crita < 3 && critd < 0.03/d && critp<0.01*maxRatio)
        break;
			iters+=1;
    }
    return crita < 3 && critd < 0.03/d;//true;
  }

	bool ok=true;

  ok=solve();
/*
  a=239.89983799131;
	d=131.208741788796;
	phase=0.471433072970658;
*/

  write("a: ", a);
  write("d: ", d);
  write("p: ", phase);
  write(ok ? "ok":"failed");

  return ok;

}


real SHIFT_FACTOR=3.3;

picture p;
for(int j; j<Nsigs; ++j)
{
  guide g0;
  for(pair p : signals[j])
    g0=g0--(p.x, p.y + zs[j]*SHIFT_FACTOR);
  draw(p, g0, grey+opacity(0.3));
}

int[] breezers;

for(int i=0;i<61;++i)
  breezers.push(i);

file brparams = output("breezer_params.dat");

pair getMax(pair[] arr)
{
  pair res;
  real absmax=0;
  for(pair v : arr)
  {
    if(abs(v.y) > absmax)
    {
      res = v;
      absmax = abs(v.y);
    }
  }
  return res;
}

for(int brId=0; brId < breezers.length; ++brId) {

	write("----------->", breezers[brId]);

  bool ok=findBreezerFor(signals[breezers[brId]], zs[breezers[brId]], false);

	if(!ok)
		continue;
/*
  a=316.08726679256;
  d=90.553106839889;
  phase=0.0286002749146886;
*/
  write("---a=",a);
  write("---d=",d);
  write("---phase=",phase);
	pair[] maxesA=getControlPointsAnalytic(phase,a,d);

	guide ga;

	pair[] maxes=getControlPoints(signals[breezers[brId]]);
	
	for(int i=0;i<maxes.length; ++i)
		maxes[i] = (maxes[i].x, maxes[i].y / beta);
	
	write("==========");
	write(maxesA);
	write(maxes);
	write("==========");	
	
  
  pair maxA = getMax(maxesA);
  pair maxsig = getMax(maxes);
  
  real x0 = -maxA.x+maxsig.x;
	real factor = maxA.y * maxsig.y < 0 ? -1 : 1;//maxes[1].y/maxesA[1].y;

	write(brparams, string(zs[brId])+" "+string(a)+" "+string(d)+" "+string(phase)+" "+string(x0)+" "+string(factor), endl);

	for(real x=x0-0.1; x<x0+0.1; x+=(signals[breezers[brId]][signals[breezers[brId]].length-1].x-signals[breezers[brId]][0].x)*0.001)
		ga=ga--(x,breezerPhase(x-x0, phase, a, d) * beta*factor + zs[breezers[brId]]*SHIFT_FACTOR);

	/*
	write(maxesA);
	for(pair pt : maxesA)
		dot(p,(pt.x+x0, pt.y*factor));
	*/

	//for(pair pt : maxes)
	//	dot(p,(pt.x, pt.y));
	//write(maxes);

	guide g0;
	for(pair p : signals[breezers[brId]])
		g0=g0--(p.x, p.y + zs[breezers[brId]]*SHIFT_FACTOR);
	
	draw(p, g0, blue);

	draw(p, ga, red);

}


guide g0;
for(pair p : signals[0])
	g0=g0--(p.x, p.y + zs[0]*SHIFT_FACTOR);
draw(p, g0, 0.5*green);


xaxis(p, "$\eta$", BottomTop, LeftTicks(0.1,0.01));
yaxis(p, "$u$", LeftRight, RightTicks);
size(p, 27cm, 5.5cm, point(p,SW), point(p,NE));
add(p.fit());

//write("zs:", zs);

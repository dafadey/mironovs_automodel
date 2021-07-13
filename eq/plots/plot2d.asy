import graph;
import palette;

currentpen = fontsize(8pt);

struct params
{
  real dz;
  int Nx;
  real Lx;
  real a;
  int skipCount;    

  real q;
  real Amp;
  real z0;
  
  void operator init(real q_, real Amp_, real z0_)
  {
    q=q_;
    Amp=Amp_;
    z0=z0_;
    dz=0.000003125;
    Nx=4096;
    Lx=1.0;
    a=0.01;
    skipCount=256;
  }

  string name()
  {
    string qs=format("%.6g", q);
    string Amps=format("%.6g", Amp);
    string z0s=format("%.6g", z0);
    
    string name="Nx="+string(Nx)+"dz="+format("%.6g",dz)+ "Lx="+format("%.6g",Lx)+"q="+qs+"a="+format("%.6g", a)+"Amp="+Amps+"z0="+z0s+".dat";
    
    return name;
  }
  
};

void addplot(picture p, string yname, string xname, pair pos, string name, pair sz=(8cm, 6cm))
{
  xaxis(p, xname, BottomTop, LeftTicks);
  yaxis(p, yname, LeftRight, RightTicks);
  p = shift(-point(p, SW)) * p;
  size(p, sz.x, sz.y, point(p, SW), point(p, NE));
  add(p.fit(), pos);
  label(name, pos+(sz.x/2, sz.y),N);
}

pen[] myCB(real min, real max)
{
  write("myCB (min, max) = ",(min,max));
  int maxn=1024;
  real absmax=max(abs(min), abs(max));
  real dv = absmax*2/maxn;
  int n = floor((max-min)/dv);
  pen[] res = new pen[n];
  write("myCB n = ",n);
  for(int i=0; i < n ; ++i)
  {
    real v = (min + i/(n-1)*(max-min))/absmax;
    res[i] = (v>0 ? (abs(v))^0.4*red : 0*red) + (v<0 ? (abs(v))^0.4*blue : 0*blue) + (abs(v))^1.5*green;
  }
  return res;
}

void addraw(picture pic_soliton_params, picture picsolitons, picture piccmint, picture piccm, picture pic, picture pica, picture picaw, params par) {
  file in=input("../raw_"+par.name());
  write("reading raw file");
  int Nx = par.Nx;
  real Lx = par.Lx;

  real[] raw = in;
  write("read raw file");

  int Ny = floor(raw.length/Nx);

  real DZ = par.dz * par.skipCount;
  real z0 = par.z0;
  real z1 = par.z0 + DZ * Ny;

  real rmin = 0;
  real rmax = 0;

  for(int i=0;i<Nx; ++i)
  {
    rmin = min(rmin, raw[i]);
    rmax = max(rmax, raw[i]);
  }

  real factor = 1/(rmax-rmin) * abs(z1-z0) * 0.06;


  write("reading momenta file");
  file in_momenta=input("../momenta_"+par.name());
  real[] raw_momenta = in_momenta;
  int momenta_items_per_row=floor(raw_momenta[0]);
  write("read momenta file");
  int Ny_momenta = floor((raw_momenta.length-1)/momenta_items_per_row);
  guide gm;
  picture pic_tmp;

  pair[] wshade;

  int jmin = floor(Nx/4);
  int jmax=floor(Nx/4*3);
  write("****Ny=", Ny);
  {
  file fout=output("signal.dat");
  int Nsignals=61;
  write(fout, string(Nsignals), endl);
   for(int sigId=0;sigId<Nsignals;++sigId)
   {
	  write(fout, string(Nx)+" "+string(z0+(floor(Ny/5)+32*sigId)/Ny*(z1-z0)), endl);
	  for(int j=0; j<Nx; ++j)
		write(fout,string(j/Nx * Lx)+" "+string(raw[j+(floor(Ny/5)+32*sigId)*Nx]), endl);
	}  
  }

  write("starting to plot simulated data");
  int lastNy=0;
  for(int i=0; i<Ny; i+=48*2)
  {
    if(i % floor(Ny/100) == 0)
      write(string(floor(i/Ny*100))+"%");
    real y = z0 + i/Ny * (z1-z0);
    real sum=0;
    for(int j=0; j<Nx; ++j)
      sum+=abs(raw[j+i*Nx]);
    sum/=Nx;
    real ssum=0;
    for(int j=Nx-113; j<Nx; ++j)
      ssum+=abs(raw[j+i*Nx]);
    ssum/=113;
    if(ssum<0.01*sum)
    {
      lastNy=i;
      guide g;
      for(int j=jmin; j<jmax; ++j)
        g=g--(j/Nx * Lx, y + factor * raw[j+i*Nx]);
	  draw(pic_tmp, (jmin/Nx * Lx,y)--(jmax/Nx * Lx,y), grey);
      draw(pic_tmp, g);
      if(i<Ny_momenta) {
        gm=gm--(raw_momenta[i*momenta_items_per_row+3+1],y);
        wshade.push((raw_momenta[i*momenta_items_per_row+3+1]-raw_momenta[i*momenta_items_per_row+4+1],y));
        wshade.push((raw_momenta[i*momenta_items_per_row+3+1]+raw_momenta[i*momenta_items_per_row+4+1],y));
      }
    }
  }
  
  real[][] img = new real[jmax-jmin][lastNy];
  real[][] img2 = new real[jmax-jmin][lastNy];

  for(int i=0; i<lastNy; ++i)
  {
    for(int j=0;j<jmax-jmin;++j)
      img2[j][i]=0;
  }
  
  for(int i=0; i<lastNy; ++i)
  {
    for(int j=jmin;j<jmax;++j)
      img[j-jmin][i]=raw[j+i*Nx];
  }

  real absmax = max(abs(min(img)),abs(max(img)));
  for(int i=0; i<lastNy; ++i)
  {
    real absmaxline=0;
    for(int j=0;j<jmax-jmin;++j)
      absmaxline = max(abs(img[j][i]), absmaxline);
    
    for(int j=1;j<jmax-jmin-1;++j)
    {
      if(abs(img[j][i])>0.1*absmax)
      {
        real cr = (img[j-1][i]+img[j+1][i]-2*img[j][i])/img[j][i];
        if(cr > 0)
          img2[j][i]=1;
        if(cr < 0)
          img2[j][i]=-1;
      }
    }
  }

  pair[][] soliton_traces;
  for(int i=0; i<lastNy; ++i)
  {
    pair[] solitons;
    real absmaxv=0;
    int j_at_max=0;
    bool in=false;
    bool soliton=false;
    for(int j=1;j<jmax-jmin-1;++j)
    {
      if(img2[j][i]!=0 && in==false)
      {
        soliton=false;
        absmaxv=0;
        in=true;
      }
      if(img2[j][i]==0 && in==true)
      {
        if(soliton && abs(-2*img[j_at_max][i]+img[j_at_max+1][i]+img[j_at_max-1][i])>0.0001*absmax) {
          //y=y0+bx+ax^2
          //y(+1) = y0+b+a
          //y(-1) = y0-b+a
          // a = (y(+1)+y(-1)-2y0)/2
          real a=(-2*img[j_at_max][i]+img[j_at_max+1][i]+img[j_at_max-1][i])/2;
          // b = (y(+1)-y(-1))/2
          real b=(img[j_at_max+1][i]-img[j_at_max-1][i])/2;
          // xmax:   y`=b+2ax ->   x = -b/2a
          real xmax = -b/2/a;
          // ymax=y0+bxmax+axmax^2
          real ymax = img[j_at_max][i] + b*xmax+a*xmax^2;
          solitons.push((j_at_max+xmax+jmin, ymax));
        }
        in=false;
        soliton=false;
      }
      if(in)
      {
        if(img2[j][i]==1)
          soliton=true;
        if(abs(img[j][i])>absmaxv)
        {
          absmaxv = abs(img[j][i]);
          j_at_max = j;
        }
      }
      
    }
    soliton_traces.push(solitons);
  }
  
  image(piccm, img, (jmin*par.Lx/par.Nx,par.z0), (jmax*par.Lx/par.Nx,par.z0+lastNy*DZ), myCB(min(img), max(img)));


  real[][] imgint = new real[jmax-jmin][lastNy];
  
  for(int i=0; i<lastNy; ++i)
  {
    for(int j=0;j<jmax-jmin;++j)
      imgint[j][i]=0;
  }
  
  for(int i=0; i<lastNy; ++i)
  {
    for(int j=jmax-jmin-2;j>0;--j)
      imgint[j][i]=imgint[j+1][i]+img[j][i];
  }
  
  image(piccmint, imgint, (jmin*par.Lx/par.Nx,par.z0), (jmax*par.Lx/par.Nx,par.z0+lastNy*DZ), myCB(min(imgint), max(imgint)));


  triple[][] traces;
  for(int i=0; i<lastNy; ++i)
  {
    for(pair v : soliton_traces[i])
    {
      bool found=false;
      for(triple[] trace : traces)
      {
        if(abs(trace[trace.length-1].x-v.x)<7 && i-trace[trace.length-1].y <= 3)
        {
          trace.push((v.x,i,v.y));
          found = true;
        }
      }
//      write("starting new trace " , (v,i));
      if(!found)
        traces.push(new triple[]{(v.x,i,v.y)});
    }
  }

  real maxW=0;
  real maxA=0;
  
  for(triple[] trace : traces)
  {
    if(trace.length<13)
      continue;
    write("drawing trace of length ", trace.length);
    guide gt;
    guide gtt;
    guide ga;
    guide gaa;
    guide gaw;
    pen tp = rand()/randMax*red+rand()/randMax*blue+rand()/randMax*green;

    for(int i=1; i<trace.length-1; ++i) {
      triple p=trace[i];
      triple p_1=trace[i-1];
      triple p1=trace[i+1];
      real z=p.y*DZ+par.z0;
      real z_1=p_1.y*DZ+par.z0;
      real z1=p1.y*DZ+par.z0;
      real A=p.z;
      real xtrace=p.x*par.Lx/par.Nx;
      real xtrace_1=p_1.x*par.Lx/par.Nx;
      real xtrace1=p1.x*par.Lx/par.Nx;
      real v = (xtrace1-xtrace_1)/(z1-z_1);
      gt = gt -- ((xtrace, z));
      gtt = gtt -- (z, xtrace);
      ga = ga -- (z, A);
      if(v>0)
	  {
		real a = sqrt(v)*abs(z);
		maxA = max(maxA, a);   
		gaa = gaa -- (z, (A > 0 ? 1:-1) * a);
	  }
      if(v>0)
	  {
		real W = sqrt(v)*z^2*par.q;
		maxW = max(maxW, W);    
      }
    }

    draw(piccm, gt, tp);
    draw(pic, gt, tp);
    draw(picsolitons, gtt);
    draw(pic_soliton_params, ga);
    if(length(gaa) > 0)
      draw(pic_soliton_params, gaa, red);
  }

  write("maxA=", maxA);
  write("maxW=", maxW);

  for(triple[] trace : traces)
  {
    if(trace.length<13)
      continue;
    write("drawing trace of length ", trace.length);
    guide gaw;

    for(int i=1; i<trace.length-1; ++i) {
      triple p=trace[i];
      triple p_1=trace[i-1];
      triple p1=trace[i+1];
      real z=p.y*DZ+par.z0;
      real z_1=p_1.y*DZ+par.z0;
      real z1=p1.y*DZ+par.z0;
      real xtrace_1=p_1.x*par.Lx/par.Nx;
      real xtrace1=p1.x*par.Lx/par.Nx;
      real v = (xtrace1-xtrace_1)/(z1-z_1);
      real A=p.z;
      if(v>0) {
		real W = sqrt(v)*z^2*par.q;
        gaw = gaw -- (z, (A > 0 ? 1:-1) * W/maxW*maxA);
	  }
	}

    if(length(gaw) > 0)
      draw(pic_soliton_params, gaw, .5 * green);
  }



  
  guide gshade;
  for(int i=0; i<wshade.length; i+=2)
    gshade=gshade--wshade[i];
  for(int i=wshade.length-1; i>0; i-=2)
    gshade=gshade--wshade[i];
  gshade=gshade--cycle;
  fill(pic, gshade, white*0.25+opacity(0.3));
  draw(pic, gm, dashed);
  add(pic, pic_tmp);

  write("ploted simulated data");

  //analytic:
  guide gcm_sim;
  guide gcm_a;
  guide gW_sim;
  guide gWloss_sim;
  guide gWtotal_sim;
  real dx = Lx/Nx;
  real cma = raw_momenta[3+1]*raw_momenta[0+1];

  for(int i=0; i<min(Ny_momenta-13, lastNy); ++i)
  {
    if(i % floor(Ny_momenta/100) == 0)
      write(string(floor(i/Ny_momenta*100))+"%");
    real z = z0 + i/Ny * (z1-z0);
    real int1=0;
    for(int j=1;j<Nx-1;++j)
      int1 += ((raw[j+1+i*Nx]-raw[j-1+i*Nx])/2/dx)^2*dx;
    real int2=0;
    for(int j=1;j<Nx-1;++j)
      int2 += raw[j+i*Nx]^4*dx;
    cma += (-3*par.q*int1 + 3*int2/z^2)*DZ;
    gcm_a=gcm_a--(z,cma/raw_momenta[0+1] - raw_momenta[3+1]);
    gcm_sim=gcm_sim--(z,raw_momenta[momenta_items_per_row*i+3+1] - raw_momenta[3+1]);
    gW_sim=gW_sim--(z,raw_momenta[momenta_items_per_row*i+1]);
    if(momenta_items_per_row >= 6)
    {
      gWloss_sim=gWloss_sim--(z,raw_momenta[momenta_items_per_row*i+5+1]);
      gWtotal_sim=gWtotal_sim--(z,raw_momenta[momenta_items_per_row*i+1]+raw_momenta[momenta_items_per_row*i+5+1]);
    }
  }
  write("ploted analytic data");

  write("gcm_a.length()=",length(gcm_a));
  write("gcm_sim.length()=",length(gcm_sim));
  write("gW_sim.length()=",length(gW_sim));
  write("gWloss_sim.length()=",length(gWloss_sim));
  write("gWtotal_sim.length()=",length(gWtotal_sim));
  
  draw(pica,gcm_a,dashed);
  draw(pica,gcm_sim);
  draw(picaw,gW_sim,0.6*green);
  if(length(gWloss_sim) > 0)
    draw(picaw,gWloss_sim,0.7*red);
  if(length(gWtotal_sim) > 0)
    draw(picaw,gWtotal_sim, gray);
  ylimits(picaw,0);
}



void add(params par, pair pos)
{
  picture picsolitonparams;
  picture picsolitons;
  picture piccm;
  picture piccmint;
  picture pic;
  picture pica;
  picture picaw;
  addraw(picsolitonparams, picsolitons, piccmint, piccm, pic, pica, picaw, par);
  addplot(piccm, "$z$", "$\eta$", pos, "q="+format("%.3g",par.q)+" Amp="+format("%.3g",par.Amp)+" Nx="+string(par.Nx), (6.7cm,6.7cm));
  addplot(piccmint, "$z$", "$\eta$", pos-(0, 7cm), "q="+format("%.3g",par.q)+" Amp="+format("%.3g",par.Amp)+" Nx="+string(par.Nx), (6.7cm,6.7cm));

  addplot(picsolitonparams, "$A, W$", "$z$", pos-(24cm, 0), "Soliton parameters", (6.5cm,6.7cm));

  addplot(picsolitons, "$\eta$", "$z$", pos-(16cm, 0), "Soliton trajectories", (6.5cm,6.7cm));

  addplot(pic, "$z$", "$\eta$", pos-(8cm,0), "q="+format("%.3g",par.q)+" Amp="+format("%.3g",par.Amp)+" Nx="+string(par.Nx), (7.7cm,6.7cm));
  addplot(pica, "$\eta$", "$z$", pos+(8.1cm, 0), "Mass Ceneter (analytic is dashed)", (5cm,6.7cm));
  addplot(picaw, "W", "$z$", pos+(14.5cm, 0), "energy",(3.5cm,6.7cm));

}


//add(params(/*q*/ 0.0000004, /*Amp*/ 0.2, /*z0*/ 1), (0,0cm));
//add(params(/*q*/ 0.0000004, /*Amp*/ 0.4, /*z0*/ 1), (0,7.5cm));
//add(params(/*q*/ 0.0000004, /*Amp*/ 0.6, /*z0*/ 1), (0,15cm));
//add(params(/*q*/ 0.0000004, /*Amp*/ 0.8, /*z0*/ 1), (0,22.5cm));

//add(params(/*q*/ 0.0000002, /*Amp*/ 0.5, /*z0*/ 1), (0,0cm));
//add(params(/*q*/ 0.0000004, /*Amp*/ 0.5, /*z0*/ 1), (0,10cm));
//add(params(/*q*/ 0.0000008, /*Amp*/ 0.5, /*z0*/ 1), (0,20cm));

//add(params(/*q*/ 0.00000083, /*Amp*/ 0.5, /*z0*/ 1), (0,0cm));
//add(params(/*q*/ 0.000000833, /*Amp*/ 0.5, /*z0*/ 1), (0,7.5cm));


//add(params(/*q*/ 0.0000008, /*Amp*/ 0.2, /*z0*/ -2), (0,0cm));
//add(params(/*q*/ 0.0000008, /*Amp*/ 0.4, /*z0*/ -2), (0,7.5cm));
//add(params(/*q*/ 0.0000008, /*Amp*/ 0.8, /*z0*/ -2), (0,15cm));

//add(params(/*q*/ 0.0000002, /*Amp*/ 0.2, /*z0*/ -2), (0,0cm));
//add(params(/*q*/ 0.0000004, /*Amp*/ 0.2, /*z0*/ -2), (0,7.5cm));
//add(params(/*q*/ 0.0000008, /*Amp*/ 0.2, /*z0*/ -2), (0,15cm));


//add(params(/*q*/ 0.000000813, /*Amp*/ 0.5, /*z0*/ 1), (0,0));

//add(params(/*q*/ 0.000000823, /*Amp*/ 0.17, /*z0*/ 1), (0,0));
//add(params(/*q*/ 0.000000823, /*Amp*/ 0.2, /*z0*/ 1), (0,7.5cm));
//add(params(/*q*/ 0.000000823, /*Amp*/ 0.25, /*z0*/ 1), (0,15cm));

add(params(/*q*/ 0.0000004, /*Amp*/ 0.6, /*z0*/ 1), (0,0));
//add(params(/*q*/ 0.000000888, /*Amp*/ 0.5, /*z0*/ 1), (0,0));


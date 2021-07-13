#include <iostream>
#include <simpledraw.h>
#include <fft_real_3.h>
#include <vector>
#include <math.h>
#include <fstream>
#include <sstream>

enum eTYPE{EVEN, ODD};

template <typename T>
struct model {
  T Wloss{};
  int Nx;
  T Lx;
  T q;
  T dz;
  T z;
  std::vector<T> vdata;
  std::vector<T> kdata;
  std::vector<T> k3data;
  std::vector<T> k3data_corrector;
  std::vector<T> tmpdata;
  void init(T a, T A, eTYPE type = eTYPE::EVEN);
  T step();
  T cm() const;
  T width() const;
  T S() const;
  T absS() const;
  T W() const;
  bool zconst{};
  std::vector<T> dvdz;
  std::vector<T> dvdx;
  T xcp{0.6};
};

template <typename T>
void model<T>::init(T a, T A, eTYPE type) {
  if(vdata.size() != Nx)
    vdata.resize(Nx, 0);
  if(kdata.size() != Nx)
    kdata.resize(Nx);
  if(k3data.size() != Nx)
    k3data.resize(Nx);
  if(k3data_corrector.size() != Nx)
    k3data_corrector.resize(Nx);
  if(tmpdata.size() != Nx)
    tmpdata.resize(Nx);
  if(dvdz.size() != Nx)
    dvdz.resize(Nx);
  if(dvdx.size() != Nx)
    dvdx.resize(Nx);
    
  for(int i=0; i<Nx; i++) {
    T x = (T) i * Lx / (T) Nx;
    T xl  = (x - Lx * xcp) / a;
    vdata[i] = A * xl * exp( - xl * xl);
  }
  if(type == eTYPE::ODD)
  {
    qqFFT_freal<T>(Nx, vdata.data(), tmpdata.data());

    for(int i=0; i<Nx; i += 2) {
      int ie = i;
      int io = i + 1;
      const T vc = vdata[ie];
      const T vs = vdata[io];
      const T k = (T) i * M_PI / Lx;

      vdata[ie] = k * vs;
      vdata[io] = -k * vc;
    }    

    qqFFT_freal_1<T>(Nx, vdata.data(), tmpdata.data());

    T vmax=0;
    for(const auto& v : vdata)
      vmax = std::max(std::abs(v), vmax);
    for(auto& v : vdata)
      v = v / vmax * A;
  }
}

template <typename T>
T model<T>::step() {
  for(int i=0; i<Nx; i++)
    k3data[i] = vdata[i] * vdata[i] * vdata[i];
  qqFFT_freal<T>(Nx, vdata.data(), tmpdata.data());
  qqFFT_freal<T>(Nx, k3data.data(), tmpdata.data());
  T sf = vdata[0];
  
  int ii=0;
  
  static std::vector<T> cos3s;
  static std::vector<T> sin3s;
  
  if(cos3s.size() == 0)
  {
    cos3s.resize(Nx/2);
    sin3s.resize(Nx/2);
    cos3s[0] = 1;
    const T k = (T) (Nx) * M_PI / Lx;
    const T k3 = k * k * k;
    sin3s[0] = sin(-q * k3 * dz);
    int ii=0;
    for(int i=2; i<Nx; i += 2, ii += 1) {
      const T k = (T) i * M_PI / Lx;
      const T k3 = k * k * k;
      cos3s[ii] = cos(-q * k3 * dz);
      sin3s[ii] = sin(-q * k3 * dz);
    }    
  }

  for(int i=0, ii=0; i<Nx; i += 2, ii += 1) {
    int ie = i;
    int io = i + 1;
    const T vc = vdata[ie];
    const T vs = vdata[io];
    const T k = (T) i * M_PI / Lx;
    const T k3 = k * k * k;
//    const T z_2 = z>1.8 ? 0 : 1./(z * z);
    const T z_2 = zconst ? 1. : 1./(z * z);

    vdata[ie] = cos3s[ii] * vc - sin3s[ii] * vs + (-2.) * z_2 * k * k3data[io] * dz;
    vdata[io] = cos3s[ii] * vs + sin3s[ii] * vc - (-2.) * z_2 * k * k3data[ie] * dz;
    
    dvdz[ie] = vdata[ie] - vc;
    dvdz[io] = vdata[io] - vs;
    
    dvdx[ie] = -vs * k;
    dvdx[io] = vc * k;
    //note first two members are due to exact linear solution of dv/dz=-ik^3 v -> v1=v0*exp(-ik^3 * dz)
    // second member is added with simple euller methind with fourier approximation of derivative along x
  }

  qqFFT_freal_1<T>(Nx, vdata.data(), tmpdata.data());
  //find corrector
  for(int i=0; i<Nx; i++)
    k3data_corrector[i] = vdata[i] * vdata[i] * vdata[i];
  qqFFT_freal<T>(Nx, k3data_corrector.data(), tmpdata.data());
  qqFFT_freal<T>(Nx, vdata.data(), tmpdata.data());
  
  for(int i=0,ii=0; i<Nx; i += 2, ii += 1) {
    int ie = i;
    int io = i + 1;
    const T vc = vdata[ie];
    const T vs = vdata[io];
    const T k = (T) i * M_PI / Lx;
    const T k3 = k * k * k;

//    const T z_2 = z>1.8 ? 0 : 1./(z * z);
//    const T z_2n = z>1.8 ? 0 : 1./((z+dz) * (z+dz));
    const T z_2 = zconst ? 1. : 1./(z * z);
    const T z_2n = zconst ? 1. : 1./((z+dz) * (z+dz));
    
    vdata[ie] += (-2.) * k * (-k3data[io] * z_2  + k3data_corrector[io] * z_2n) * 0.5 * dz;
    vdata[io] += - (-2.) * k * (-k3data[ie] * z_2 + k3data_corrector[ie] * z_2n) * 0.5 * dz;

    //note first two members are due to exact linear solution of dv/dz=-ik^3 v -> v1=v0*exp(-ik^3 * dz)
    // second member is added with simple euller methind with fourier approximation of derivative along x
  }
  
  qqFFT_freal_1<T>(Nx, vdata.data(), tmpdata.data());


  qqFFT_freal_1<T>(Nx, dvdz.data(), tmpdata.data());
  qqFFT_freal_1<T>(Nx, dvdx.data(), tmpdata.data());

  
  static std::vector<T> mask;
  if(mask.size() == 0)
  {
    mask.resize(Nx);
    for(int i=0;i<Nx;i++)
    {
      T x = (T) i / (T) Nx;
      mask[i] = 0.5 - 0.5 * tanh((x-0.05)*(0.95-x)*30.0);
    }
  } 
  
  T W0 = W();
  for(int i=0;i<Nx;i++)
  {
    T x = (T) i / (T) Nx;
    vdata[i] *= 1.0 - dz * /*39.065*/ 100 * mask[i];
  }
  
  T W1 = W();
  
  Wloss += W0-W1;
  
  z += dz;  

  // dv/dz + q*d3 v/dx^3 - 2 / z^2 dv^3/dx = 0
  // dp/dz - i k^3 q * p - 2 / z^2 i k F(v^3) = 0
  return sf;
}

template <typename T>
T model<T>::S() const {
  T s{0.};
  T dx = Lx / (T) Nx;
  for(int i=0; i < Nx; i++)
    s += vdata[i] * dx;
  return s;
}

template <typename T>
T model<T>::absS() const {
  T s{0.};
  T dx = Lx / (T) Nx;
  for(int i=0; i < Nx; i++)
    s += fabs(vdata[i]) * dx;
  return s;
}

template <typename T>
T model<T>::W() const {
  T s{0.};
  T dx = Lx / (T) Nx;
  for(int i=0; i < Nx; i++)
    s += vdata[i] * vdata[i] * dx;
  return s;
}

template <typename T>
T model<T>::cm() const {
  T norma{0.};
  T cm{0.};
  T dx = Lx / (T) Nx;
  for(int i=0; i < Nx; i++) {
    T x = (T) i / (T) Nx * Lx;
    cm += vdata[i] * vdata[i] * x * dx;
    norma += vdata[i] * vdata[i] * dx;
  }
  return cm / norma;
}

template <typename T>
T model<T>::width() const {
  T norma{0.};
  T width{0.};
  T x0 = cm();
  T dx = Lx / (T) Nx;
  for(int i=0; i < Nx; i++) {
    T x = (T) i / (T) Nx * Lx;
    width += vdata[i] * vdata[i] * (x - x0) * (x - x0) * dx;
    norma += vdata[i] * vdata[i] * dx;
  }
  return sqrt(width / norma);
}

int main(int argc, char* argv[])
{
  #define FL_DBL double
  int Nx = 512;
  FL_DBL Lx=0.5;
  FL_DBL dz=0.000003;
  FL_DBL z0=1.;
  FL_DBL z1=3.;
  FL_DBL q=0.0001/4;
  FL_DBL a=0.01;
  FL_DBL Amp=5.64;
  int skipCount = 8;
  int type = 0;
  int store = 1;
  int draw=1;
  int zconst=0;
  FL_DBL xcp=0.6;
  std::string prefix{};
  #define PARSE(ARG) if(name == #ARG) { sval >> ARG; continue;}
  for(int i=1;i<argc;i++)
	{
		std::string inp = std::string(argv[i]);
		size_t pos = inp.find("=");
		if(pos == std::string::npos)
			printf("you specified parameter wrong way use <name>=<value> format. NOTE: no \'-\' and spaces\n");
		else
		{
			std::string name = inp.substr(0,pos);
			std::stringstream sval;
			sval << inp.substr(pos+1,std::string::npos);
			printf("parameter[%d] has name %s and value %s\n",i-1,name.c_str(), sval.str().c_str());
			PARSE(Nx);
			PARSE(Lx);
			PARSE(dz);
			PARSE(z0);
			PARSE(z1);
			PARSE(q);
			PARSE(a);
			PARSE(Amp);
			PARSE(skipCount);
			PARSE(type);
			PARSE(store);
			PARSE(draw);
			PARSE(zconst);
			PARSE(xcp);
			PARSE(prefix);
    }
  }

  std::cout << "Hi!\n";
  model<FL_DBL> v;
  v.Nx=Nx;
  v.Lx=Lx;
  v.dz=dz;
  v.q=q;
  v.z=z0;
  v.xcp = xcp;
  v.init(a, Amp, type == 1 ? eTYPE::ODD : eTYPE::EVEN);
  v.zconst = zconst;
  FL_DBL W0 = v.W();

  std::stringstream basenamess;
  basenamess << "Nx=" << v.Nx << "dz=" << v.dz << "Lx=" << v.Lx << "q=" << v.q << "a=" << a << "Amp=" << Amp << "z0=" << z0 << ".dat";
  std::string basename=basenamess.str();

  std::string outname = prefix + "momenta_" + basename;
 
  std::string rawname = prefix + "raw_" + basename;

  std::ofstream output;
  if(store)
    output.open(outname.c_str());
  
  output << 6 << '\n'; // number of items per row
  
  int steps=floor((z1-z0)/dz);//128*8*4;

  std::cout << "steps=" << steps << '\n';
  std::vector<FL_DBL> display(v.Nx * ceil((double)steps / (double)skipCount) , 0.);
  if(draw)
    fadey_init(v.Nx, display.size() / v.Nx, 1);

  for(int i = 0; i < steps; i++)
  {
    auto sf = v.step();
    if(i % skipCount == 0)
    {
      std::cout << "z=" << v.z << ", W/W0=" << v.W()/W0 << ", abs=" << v.absS() << ", sf=" << sf << ", s=" << v.S() << ", cm=" << v.cm() << ", w=" << v.width() << '\n';
      if(store)
        output << v.W() << "\t" << sf << "\t" << v.S() << "\t" << v.cm() << "\t" << v.width() << "\t" << v.Wloss << '\n';
      for(int j = 0; j < v.Nx; j++)
        display[floor(i / skipCount) * v.Nx + j] = v.vdata[j];
      if(draw)
        fadey_draw(display.data(), v.Nx, display.size() / v.Nx, 0);
    
      for(int j=v.Nx/4; j<v.Nx*3/4; j++)
      {
        if(std::abs(v.dvdz[j]) > std::abs(v.dvdx[j]))
          std::cout << j << " " << v.dvdz[j] << " : " << v.dvdx[j] << '\n';
      }

    }
    bool nan=false;
    for(int j = 0; j < v.Nx; j++)
      nan = v.vdata[j] == v.vdata[j] ? nan : true;

    if(nan)
      break;
  }

  if(store) {
    output.close();
    
    std::cout << "saving raw\n";
    
    std::ofstream rawoutput(rawname.c_str());
    for(int i=0; i < display.size() / v.Nx; i++) {
      for(int j=0; j<v.Nx; j++)
        rawoutput << display[j+i*v.Nx] << '\n';
    }
    rawoutput.close();
  }
  
  std::cout << "done\n";
 
  /*
  char a;
  std::cin >> a;
  */
  if(draw)
    fadey_close();
  return 0;
}

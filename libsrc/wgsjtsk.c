// Zjednodušený převod polohových souřadnic mezi systémem WGS-84 a S-JTSK

// Geodetické referenční systémy v České republice : Vývoj od klasických ke geocentrickým souřadnicovým systémům.
// Odbor. red. Jan Kostelecký, Drahomír Dušátko. Praha : VZÚ, 1998.
// Praha (Zdiby): VÚGTK, 1998, Edice: Edice VÚGTK Ročník 44


void wgs2jtsk (double phi, double lambda, double *x, double *y ) {
// Převod WGS-84 -> S-JTSK
// Převod je určen pro souřadnice v ČR

const double p0 = 50;
const double l0 = 15;

const double y0 = 703000;
const double x0 = 1058000;

const double Ay =  1.180672981e+01; // 1
const double By = -1.431119075e+04; // p
const double Cy = -7.109369068e+04; // l
const double Dy =  4.527213114e-02; // p*p
const double Ey =  1.469297520e+03; // p*l
const double Fy = -6.216573827e+01; // l*l
const double Gy =  1.746024222e+00; // p*p*p
const double Hy =  1.482366057e+00; // p*p*l
const double Iy = -1.646574057e+00; // p*l*l
const double Jy =  1.930950004e+00; // l*l*l

const double Ax =  1.471808238e+02;  // 1
const double Bx = -1.102950611e+05;  // p
const double Cx =  9.224512054e+03;  // l
const double Dx = -1.335425822e+01;  // p*p
const double Ex = -1.928902631e+02;  // p*l
const double Fx = -4.735502716e+02;  // l*l
const double Gx = -4.564660084e+00;  // p*p*p
const double Hx = -4.355296392e+00;  // p*p*l
const double Ix =  8.911019558e+00;  // p*l*l
const double Jx =  3.614170182e-01;  // l*l*l

double p,l;
p=phi-p0;
l=lambda-l0;

*x = x0 + Ax + Bx * p + Cx * l + Dx * p*p + Ex * p*l + Fx * l*l + Gx * p*p*p + Hx * p*p*l + Ix * p*l*l + Jx * l*l*l;
*y = y0 + Ay + By * p + Cy * l + Dy * p*p + Ey * p*l + Fy * l*l + Gy * p*p*p + Hy * p*p*l + Iy * p*l*l + Jy * l*l*l;

}


void jtsk2wgs (double X, double Y, double *phi, double *lambda ) {
// Převod S-JTSK -> WGS-84
// Převod je určen pro souřadnice v ČR

const double p0 = 50;
const double l0 = 15;

const double y0 = 703000;
const double x0 = 1058000;

const double Ap =  1.325132993e-03;  // 1
const double Bp = -8.916429099e-06;  // x
const double Cp = -1.156917384e-06;  // y
const double Dp = -2.298750250e-14;  // x*x
const double Ep =  2.087176527e-13;  // x*y
const double Fp = -8.219794748e-13;  // y*y
const double Gp =  2.191874854e-20;  // x*x*x
const double Hp =  5.305545189e-21;  // x*x*y
const double Ip =  1.760134043e-19;  // x*y*y
const double Jp =  6.270628603e-21;  // y*y*y

const double Al = -1.019442857e-04;  // 1
const double Bl =  1.794902692e-06;  // x
const double Cl = -1.383338939e-05;  // y
const double Dl = -3.294257309e-13;  // x*x
const double El =  2.506009659e-12;  // x*y
const double Fl =  3.291143794e-13;  // y*y
const double Gl =  4.567560092e-20;  // x*x*x
const double Hl = -4.843979237e-19;  // x*x*y
const double Il = -1.182561606e-19;  // x*y*y
const double Jl =  1.641107774e-19;  // y*y*y

double x,y;
x = X - x0;
y = Y - y0;

*phi    = p0 + Ap + Bp * x + Cp * y + Dp * x*x + Ep * x*y + Fp * y*y + Gp * x*x*x + Hp * x*x*y + Ip * x*y*y + Jp * y*y*y;
*lambda = l0 + Al + Bl * x + Cl * y + Dl * x*x + El * x*y + Fl * y*y + Gl * x*x*x + Hl * x*x*y + Il * x*y*y + Jl * y*y*y;

}

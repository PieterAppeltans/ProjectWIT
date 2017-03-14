#ifndef constants
#define constants

/* Global variables for every storage scenario */
const double DU_R = 2.8*pow(10, -10);
const double DU_Z = 1.1*pow(10, -9);
const double DV_R = 2.32*pow(10, -9);
const double DV_Z = 6.97*pow(10, -9);
const double rq = 0.97;

/* Global variables for applicable storage scenario */
const double ETA_U = 0.208;
const double ETA_V = 0.0004;
const double T_CEL = 298.15;

const double vamb = (101300*ETA_U)/(8.314*T_CEL);
const double uamb = (101300*ETA_V)/(8.314*T_CEL);
const double Vmfv = 1.61*pow(10, -4)*exp((56700/8.314)*(1/293.15 - 1/T_CEL));
const double Vmu = 2.39*pow(10, -4)*exp((80200/8.314)*(1/293.15 - 1/T_CEL));
const double hu = 7*pow(10,-7);
const double hv = 7.5*pow(10,-7);
const double Kmu = 0.1149;
const double Kmv = 27.2438;

vector<double> Ru_simple(vector<double> u) {
  /* To calculate starting value of u for solving nonlinear system of equations */
  return (Vmu/0.4103)*u;
}

vector<double> Rv_simple(vector<double> u) {
  /* To calculate starting value of v for solving nonlinear system of equations */
  return rq*Ru_simple(u);
}

vector<double> Ru(vector<double> u, vector<double> v) {
  return element_div(Vmu*u,element_prod(scalar_vector<double>(u.size(),0.4103)+u,scalar_vector<double>(u.size(),1)+(v/27.2438)));
}

vector<double> Rv(vector<double> u, vector<double> v) {
  return 0.97*Ru(u,v)+ element_div(Vmfv*scalar_vector<double>(u.size(),1),scalar_vector<double>(u.size(),1)+(u/0.1149));
}

matrix<double> diagonalize(vector<double> v) {
  matrix<double> func(v.size(), v.size());
  for (unsigned l = 0; l < v.size(); ++l) {
    func(l, l) = v(l);
  }
  return func;
}

matrix<double> dRudu(vector<double> u, vector<double> v) {
  scalar_vector<double> one(u.size(),1.);
  vector<double> res_vec = element_div(Vmu*0.4103*one, element_prod(one+v/27.2438, element_prod(one*0.4103+u, one*0.4103+u)));
  return diagonalize(res_vec);
}

matrix<double> dRudv(vector<double> u, vector<double> v) {
  scalar_vector<double> one(u.size(),1.);
  vector<double> res_vec = element_div(-1*27.2438*Vmu*u, element_prod(0.1149*one+u, element_prod(27.2438*one+v, 27.2438*one+v)));
  return diagonalize(res_vec);
}

matrix<double> dRvdu(vector<double> u, vector<double> v) {
  scalar_vector<double> one(u.size(),1);
  vector<double> res_vec = element_div(Vmu*0.4103*one, element_prod(one+v/27.2438, element_prod(one*0.4103+u, one*0.4103+u)));
  return rq*dRudu(u, v) - diagonalize(element_div(0.1149*Vmfv*one, element_prod(0.1149*one+u, 0.1149*one+u)));
}

matrix<double> dRvdv(vector<double> u, vector<double> v) {
  return rq*dRudv(u, v);
}

#endif

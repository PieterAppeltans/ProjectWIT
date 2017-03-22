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

const double uamb = (101300*ETA_U)/(8.314*T_CEL);
const double vamb = (101300*ETA_V)/(8.314*T_CEL);
const double Vmfv = 1.61*pow(10, -4)*exp((56700/8.314)*(1/293.15 - 1/T_CEL));
const double Vmu = 2.39*pow(10, -4)*exp((80200/8.314)*(1/293.15 - 1/T_CEL));
const double hu = 7*pow(10,-7);
const double hv = 7.5*pow(10,-7);
const double Kmu = 0.4103;
const double Kmv = 27.2438;
const double Kmfu = 0.1149;


typedef Array<double,Dynamic,Dynamic> ArrayXXd;
typedef Array<double,Dynamic,1> ArrayXd;

VectorXd Ru(VectorXd& u, VectorXd& v) {
  ArrayXd u_a = u.array();
  ArrayXd v_a = v.array();
  ArrayXd res_a = (Vmu*u_a)*((Kmu+u_a)*(1+v_a/Kmv)).inverse();
  VectorXd res_vec = res_a.matrix();
  return res_vec;
}
VectorXd Rv(VectorXd& u, VectorXd& v) {
  ArrayXd u_a = u.array();
  ArrayXd v_a = v.array();
  ArrayXd res_a = rq*Ru(u,v).array()+Vmfv*(1+u_a/Kmfu).inverse();
  VectorXd res_vec = res_a.matrix();
  return res_vec;
}

MatrixXd dRudu(VectorXd& u,VectorXd& v) {
  ArrayXd u_a = u.array();
  ArrayXd v_a = v.array();
  ArrayXd res_a = (Vmu*0.4103)*((1+v_a/27.2438)*(0.4103+u_a)*(0.4103+u_a)).inverse();
  VectorXd res_vec = res_a.matrix();
  DiagonalMatrix<double, Dynamic> m = res_vec.asDiagonal();
  return m;
}

MatrixXd dRudv(VectorXd& u,VectorXd& v) {
  ArrayXd u_a = u.array();
  ArrayXd v_a = v.array();
  ArrayXd res_a = (-27.2438*Vmu*u_a)*((0.1149+u_a)*(27.2438+v_a).pow(2)).inverse();
  VectorXd res_vec = res_a.matrix();
  DiagonalMatrix<double, Dynamic> m = res_vec.asDiagonal();
  return m;
}
// Na te kijkens
MatrixXd dRvdu(VectorXd& u,VectorXd& v) {
  ArrayXd u_a = u.array();
  ArrayXd v_a = v.array();
  ArrayXd res_a = rq*(Vmu*0.4103)*((1+v_a/27.2438)*(0.4103+u_a)*(0.4103+u_a)).inverse()-Vmfv*((Kmfu+u_a)*(Kmfu+u_a)).inverse();
  VectorXd res_vec = res_a.matrix();
  DiagonalMatrix<double, Dynamic> m = res_vec.asDiagonal();
  return m;
}

MatrixXd dRvdv(VectorXd& u,VectorXd& v) {
  ArrayXd u_a = u.array();
  ArrayXd v_a = v.array();
  ArrayXd res_a =  rq*(-27.2438*Vmu*u_a)*((0.1149+u_a)*(27.2438+v_a).pow(2)).inverse();
  VectorXd res_vec = res_a.matrix();
  DiagonalMatrix<double, Dynamic> m = res_vec.asDiagonal();
  return m;
}

#endif

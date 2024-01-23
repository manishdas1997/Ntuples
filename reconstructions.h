#include <cmath>
#include <array>
#include <tuple>
#include <TTree.h>
#include <TMath.h>
#include <TMatrixD.h>
#include <TVectorD.h>
#include <TRotation.h>

const double kSpeedOfLight = 29.9792458; // [cm/ns]

double d(double x, double y, double x0, double y0) {
  return sqrt(pow(x-x0, 2.) + pow(y-y0, 2.));
}

double dFx(double x, double y, double xa, double ya, double xb, double yb) {
  double ret = 0;
  ret += (x-xa) / d(x, y, xa, ya);
  ret -= (x-xb) / d(x, y, xb, yb);
  return ret;
}

double dFy(double x, double y, double xa, double ya, double xb, double yb) {
  double ret = 0;
  ret += (y-ya) / d(x, y, xa, ya);
  ret -= (y-yb) / d(x, y, xb, yb);
  return ret;
}

double F(double x, double y, double xa, double ya, double xb, double yb) {
  double ret = 0;
  ret += d(x, y, xa, ya);
  ret -= d(x, y, xb, yb);
  return ret;
}

std::tuple<int, double, TVector3> hyperbolicReco(const std::array<TVector3, 3>& hits,
                                           const std::array<double, 3>& times)
{
  int error = 0;
  
  // find the decay plane
  TVector3 normal = ((hits[1]-hits[0]).Cross(hits[2]-hits[0] )).Unit();
  
  // prepare transformation to the decay plane
  TVector3 z(0.,0.,1.);
  TVector3 rotAxis = normal.Cross( z );
  double angle = z.Angle( normal ); // radians
  
  TRotation rot;
  rot.Rotate(angle, rotAxis);
  
  // transform gamma hits to decay plane
  TVector3 gammas2D[3];
  for (int i = 0; i < 3; i++)
  {
    gammas2D[i] = rot * hits[i];
  }

  // for convenience, we copy the hits' positions to a matrix
  double hits2d[3][2];
  for (int i = 0; i < 3; i++)
  {
    hits2d[i][0] = gammas2D[i].X();
    hits2d[i][1] = gammas2D[i].Y();
  }

  // We start from the centre of the detector
  double x = 0.;
  double y = 0.;
  // increments used for iteration
  // (x,y) -> (x+dx, y+dy) in each subsequent interation
  double dx, dy;
  
  const double c = 29.9792458; // cm/ns

  // matrices for the least squares solution
  TMatrixD A(3, 2);
  TMatrixD b(3, 1);
  TMatrixD T(2, 3);
  TMatrixD PI(2, 2);
  TMatrixD ATB(2, 1);
  TMatrixD DX(2, 1);

  // Iterative solution of the overdetermined system
  double epsilon = 0.0001;

  int n_iter = 0;
  do
  {
    // fill the design (A) matrix
    A[0][0] = dFx(x, y, hits2d[0][0], hits2d[0][1], hits2d[1][0], hits2d[1][1]);
    A[0][1] = dFy(x, y, hits2d[0][0], hits2d[0][1], hits2d[1][0], hits2d[1][1]);
    A[1][0] = dFx(x, y, hits2d[1][0], hits2d[1][1], hits2d[2][0], hits2d[2][1]);
    A[1][1] = dFy(x, y, hits2d[1][0], hits2d[1][1], hits2d[2][0], hits2d[2][1]);
    A[2][0] = dFx(x, y, hits2d[0][0], hits2d[0][1], hits2d[2][0], hits2d[2][1]);
    A[2][1] = dFy(x, y, hits2d[0][0], hits2d[0][1], hits2d[2][0], hits2d[2][1]);

    //    A.Print();
    
    // fill the constant (b) vector
    b[0][0] = c * (times[0] - times[1]) - F(x, y, hits2d[0][0], hits2d[0][1], hits2d[1][0], hits2d[1][1]);
    b[1][0] = c * (times[1] - times[2]) - F(x, y, hits2d[1][0], hits2d[1][1], hits2d[2][0], hits2d[2][1]);
    b[2][0] = c * (times[0] - times[2]) - F(x, y, hits2d[0][0], hits2d[0][1], hits2d[2][0], hits2d[2][1]);

    //    b.Print();
    
    // calculate pseudoinverse
    T.Transpose(A);
    PI = T * A;
    PI.InvertFast(); // 2x2
    T.Transpose(A);  // recreate T
    ATB = T * b;     // 2x1
    // calculate corrections
    DX = PI * ATB;

    // update (x,y)
    x += DX[0][0];
    y += DX[1][0];

    if(++n_iter > 20) break;
    //    std::cout << DX[0][0] << " : " << DX[1][0] << std::endl;
    
  } while (fabs(DX[0][0]) > epsilon && fabs(DX[1][0]) > epsilon);

  //  std::cout << "Iteration terminated" << std::endl;

  // transform the result back to the 3D space
  TVector3 vertex = rot.Inverse() * TVector3(x, y, gammas2D[0].Z());
  double anh_time = times[0] - (hits[0]-vertex).Mag() / kSpeedOfLight;
  
  return std::make_tuple(error, anh_time, vertex);
}

std::tuple<int, double, TVector3> trilaterativeReco(const std::array<TVector3, 3>& hits,
                                                    const std::array<double, 3>& times){
  
  int error = 0;
  
  // find the decay plane
  TVector3 normal = (hits[1]-hits[0]).Cross(hits[2]-hits[1]).Unit();
  
  // prepare transformation to the decay plane
  TVector3 z(0.,0.,1.);
  TVector3 rotAxis = normal.Cross( z );
  double angle = z.Angle( normal ); // radians
  
  TRotation rot;
  rot.Rotate(angle, rotAxis);
  
  // transform gamma hits to decay plane
  TVector3 gammas2D[3];
  for(int i=0;i<3;i++){
    gammas2D[i] = rot * hits[i];
  }
  
  // solve in 2D
  int combs[][2] = {{0,1}, {1,2}, {0,2}};
  double M[3][3];
  double D[3];
  
  // fill the matrix and constants vector
  int i,j;
  for(int k=0;k<3;++k){ // k - rows
    i = combs[k][0];
    j = combs[k][1];
    M[k][0] = 2.*( gammas2D[i].X() - gammas2D[j].X() );
    M[k][1] = 2.*( gammas2D[i].Y() - gammas2D[j].Y() );
    M[k][2] = 2.*kSpeedOfLight*kSpeedOfLight*( times[j] - times[i] );       
    D[k] = pow(gammas2D[i].X(),2.)
      - pow(gammas2D[j].X(),2.)
      + pow(gammas2D[i].Y(),2.)
      - pow(gammas2D[j].Y(),2.)
      - kSpeedOfLight*kSpeedOfLight*pow(times[i],2.)
      + kSpeedOfLight*kSpeedOfLight*pow(times[j],2.);
  }
  
  // use analytical solutions: x = Ex*t+Fx, y=Ey*t+Fy
  double Ex,Ey,Fx,Fy;
  Ex = ( M[0][2]*M[1][1]-M[0][1]*M[1][2] )/( M[0][1]*M[1][0]-M[0][0]*M[1][1] );
  Fx = ( M[0][1]*D[1]-M[1][1]*D[0] )/( M[0][1]*M[1][0]-M[0][0]*M[1][1] );
  
  Ey = ( M[0][0]*M[1][2] - M[0][2]*M[1][0] )/( M[0][1]*M[1][0]-M[0][0]*M[1][1] );
  Fy = ( M[1][0]*D[0] - M[0][0]*D[1] )/( M[0][1]*M[1][0]-M[0][0]*M[1][1] );       
  
  // find t - using ready analytical solutions
  double a,b,cc,delta;

  int errFlag = 0;
  
  a = Ex*Ex + Ey*Ey - kSpeedOfLight*kSpeedOfLight;
  b = 2.*( Ex*(Fx-gammas2D[0].X()) + Ey*(Fy-gammas2D[0].Y()) + kSpeedOfLight*kSpeedOfLight*times[0] );
  cc = pow(Fx-gammas2D[0].X(), 2.) + pow(Fy-gammas2D[0].Y(), 2.) - kSpeedOfLight*kSpeedOfLight*pow(times[0], 2.);
  delta = b*b - 4.*a*cc;
  if( delta < 0. ){
    error = 1;
  }
  double sol_time[2];
  sol_time[0] = (-1.*b - sqrt(delta))/(2.*a);
  sol_time[1] = (-1.*b + sqrt(delta))/(2.*a);

  TVector3 sol_hit[2];
  for(int i = 0; i<2;++i){
    TVector3 sol2Dv( Ex*sol_time[i]+Fx, Ey*sol_time[i]+Fy, gammas2D[0].Z() );
    
    // transform the solution back to 3D
    sol_hit[i] =  rot.Inverse() * sol2Dv;
    
  }

  // check solution 2 for reasonability
  if( sol_time[1] < 0. || sol_time[1] > 50000.  ){
    error = 2;
  }
  if( sol_hit[1].Perp() > hits[0].Perp() ||
      sol_hit[1].Perp() > hits[1].Perp() ||
      sol_hit[1].Perp() > hits[2].Perp() ||
      fabs( sol_hit[1].Z() ) > 25.0 ){
    error = 4;
  }
  if( TMath::IsNaN( sol_time[1] ) ||
      TMath::IsNaN( sol_hit[1].X() ) ||
      TMath::IsNaN( sol_hit[1].Y() ) ||
      TMath::IsNaN( sol_hit[1].Z() )
      ){
    error = 8;
  }

  double anh_time = sol_time[1];
  
  return std::make_tuple(error, anh_time, sol_hit[1]);
}

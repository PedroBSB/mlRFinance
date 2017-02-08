#include <RcppArmadillo.h>
#include <RcppEigen.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppEigen)]]
#include <cmath>
using namespace Rcpp;

/********************************************************************************************************/
/***********************************         KERNEL FUNCTIONS      **************************************/
/********************************************************************************************************/
//https://github.com/accord-net/framework/blob/development/Sources/Accord.Statistics/Kernels/Gaussian%601.cs
//http://crsouza.com/2010/03/17/kernel-functions-for-machine-learning-applications/
//http://cseweb.ucsd.edu/~yoc002/arccos.html

double ANOVAKernel(Eigen::RowVectorXd x, Eigen::RowVectorXd y, Eigen::RowVectorXd parms)
{
  double sigma=parms(0);
  double d=parms(1);
  Eigen::RowVectorXd diff=(-sigma * (x-y).array().pow(2)).exp();
  double res=std::pow (diff.sum(), d) ;
  return(res);
}


double ArccosKernel(arma::vec x,arma::vec y,arma::vec parms)
{
  double k_xx_l = arma::sum(x % y);
  double k_yy_l = arma::sum(y % y);
  double k_xy_l = arma::sum(x % y);
  double theta_l;

  int L = parms.size();

  for (int l=1; l<=L; l++)
  {
    theta_l = std::acos(std::max(std::min(k_xy_l / std::sqrt(k_xx_l * k_yy_l),  1.0), -1.0));
    k_xy_l = std::pow(k_xx_l * k_yy_l, parms[l-1]/2) / M_PI * R::bessel_j(parms[l-1], theta_l);
    if (l < L)
    {
      k_xx_l  = std::pow(k_xx_l, parms[l-1]) / M_PI * R::bessel_j(parms[l-1], 0);
      k_yy_l  = std::pow(k_yy_l, parms[l-1]) / M_PI * R::bessel_j(parms[l-1], 0);
    }
  }
  return k_xy_l;
}


double BesselKernel(Eigen::RowVectorXd x,Eigen::RowVectorXd y,Eigen::RowVectorXd parms)
{
  double res;
  double order = parms(0);
  double degree= parms(1);
  double sigma = parms(2);
  double lim = 1.0/(R::gammafn(order+1.0)*std::pow(2.0,order));
  double comp = -1.0*(2* (x.cwiseProduct(y)).sum() - (x.cwiseProduct(x)).sum() - (y.cwiseProduct(y)).sum());
  double bkt =  sigma*std::sqrt(comp);
  if(bkt < 10e-5){
    res = lim;
  }
  else{
    res = R::bessel_j(bkt,order)*(std::pow(bkt,(-order)));
  }
  return res;
}


double CauchyKernel(Eigen::RowVectorXd x, Eigen::RowVectorXd y, Eigen::RowVectorXd parms)
{
  double sigma = parms(0);
  double den=1+((((x-y).array().pow(2)).sum())/(std::pow(sigma,2)));
  double res=1/den;
  return(res);
}


double ChiSquareKernel(Eigen::RowVectorXd x, Eigen::RowVectorXd y, Eigen::RowVectorXd parms)
{
  Eigen::RowVectorXd num=(x-y).array().pow(2);
  Eigen::RowVectorXd den= 0.5*(x+y);
  double res=1.0 - (num.cwiseQuotient(den)).sum();
  return(res);
}


double CircularKernel(Eigen::RowVectorXd x, Eigen::RowVectorXd y,Eigen::RowVectorXd parms)
{
  double sigma=parms(0);
  double norm = (x-y).norm();
  double res=0;
  if(norm>=sigma){
    res=(2.0/M_PI)*std::acos(-norm/sigma)-(2.0/M_PI)*(norm/sigma)*std::sqrt(1-std::pow((norm/sigma),2));
  }
  return(res);
}


double DirichletKernel(Eigen::RowVectorXd x,Eigen::RowVectorXd y,Eigen::RowVectorXd parms)
{
  double N=parms(0);//Dimension
  double delta=parms(1);
  double prod = 1;
  for (int i=0;i<x.size();i++)
  {
    double delta = x(i) - y(i);
    double num = std::sin((N + 0.5) * (delta));
    double den = 2.0 * std::sin(delta / 2.0);
    prod *= num / den;
  }

  return prod;
}


double ExponentialKernel(Eigen::RowVectorXd x,Eigen::RowVectorXd y,Eigen::RowVectorXd parms)
{
  double sigma = parms(0);
  double num = std::sqrt(((x-y).array().pow(2)).sum());
  double den=2*(std::pow(sigma,2));
  double res=(-(num/den));
  return(res);
}


double GaussianKernel(Eigen::RowVectorXd x, Eigen::RowVectorXd y, Eigen::RowVectorXd parms)
{
  double sigma = parms(0);
  double num = -1.0*std::pow((x-y).norm(),2);
  double den=2*(std::pow(sigma,2));
  double res=std::exp((num/den));
  return(res);
}


double GeneralizedHistogramIntersectionKernel(Eigen::RowVectorXd x,Eigen::RowVectorXd y, Eigen::RowVectorXd parms)
{
  double a=parms(0);
  double b=parms(1);
  double sum=0;
  for(int i=0;i<x.size();i++){
    sum=sum+std::min(std::pow(std::abs(x(i)),a),std::pow(std::abs(y(i)),b));
  }
  return(sum);
}


double GeneralizedTStudentKernel(Eigen::RowVectorXd x, Eigen::RowVectorXd y, Eigen::RowVectorXd parms)
{
  double d = parms(0);
  double res=1.0/(1.0+(std::pow(std::sqrt(((x-y).array().pow(2)).sum()),d)));
  return(res);
}


double HellingerKernel(Eigen::RowVectorXd x, Eigen::RowVectorXd y, Eigen::RowVectorXd parms)
{
  double r = 0;
  r=((x.cwiseProduct(y)).array().sqrt()).sum();
  return r;
}


//Only for positive values in X and Y
double HistogramIntersectionKernel(Eigen::RowVectorXd x,Eigen::RowVectorXd y,Eigen::RowVectorXd parms)
{
  double sum=0;
  for(int i=0;i<x.size();i++){
    sum=sum+std::min(x(i),y(i));
  }
  return(sum);
}


double HyperbolicTangentKernel(Eigen::RowVectorXd x,Eigen::RowVectorXd y,Eigen::RowVectorXd parms)
{
  double n=parms(0);
  double cc=parms(1);
  double cross = (x.cwiseProduct(y)).sum();
  double res=tanh(((1.0/n)*cross)+cc);
  return(res);
}


double InverseMultiquadraticKernel(Eigen::RowVectorXd x,Eigen::RowVectorXd y,Eigen::RowVectorXd parms)
{
  double cc=parms(0);
  double res=1.0/(std::sqrt(((x-y).array().pow(2)).sum() + (std::pow(cc,2))));
  return(res);
}


double LaplacianoKernel(Eigen::RowVectorXd x,Eigen::RowVectorXd y,Eigen::RowVectorXd parms)
{
  double sigma=parms(0);
  double num=std::sqrt(((x-y).array().pow(2)).sum());
  double den=sigma;
  double res = std::exp(-1.0*(num/den));
  return(res);
}


double LinearKernel(Eigen::RowVectorXd x,Eigen::RowVectorXd y,Eigen::RowVectorXd parms)
{
    double cc=parms(0);
    double cross = (x.cwiseProduct(y)).sum();
    double res = cross+cc;
    return(res);
}


double LogKernel(Eigen::RowVectorXd x, Eigen::RowVectorXd y, Eigen::RowVectorXd parms)
{
    double a=parms(0);
    double b=parms(1);
    double norm = (x-y).norm();
    double res=-1.0*std::log(std::pow(norm,a)+b);
    return(res);
}



double LogLinearKernel(Eigen::RowVectorXd x, Eigen::RowVectorXd y, Eigen::RowVectorXd parms)
{
    double d=parms(0);
    double res = (-1.0*std::log(std::pow(std::sqrt(((x-y).array().pow(2)).sum()),d)+1.0));
    return(res);
}


double MexicanHatKernel(Eigen::RowVectorXd x,Eigen::RowVectorXd y,Eigen::RowVectorXd parms)
{
    double a=parms(0);
    double res2 = ((1.0 -(x-y).array().pow(2)/std::pow(a,2)).cwiseProduct((-(x-y).array().pow(2)/(2.0*std::pow(a,2))).exp())).prod();
    return(res2);
}


double MorletKernel(Eigen::RowVectorXd x,Eigen::RowVectorXd y,Eigen::RowVectorXd parms)
{
    double a=parms(0);
    double res2 = ((5.0*(x-y)/a).array().cos()*(-(x-y).array().pow(2)/(2.0*std::pow(a,2))).exp()).prod();
    return(res2);
}


double MultiquadraticKernel(Eigen::RowVectorXd x, Eigen::RowVectorXd y, Eigen::RowVectorXd parms)
{
    double cc=parms(0);
    double res = std::sqrt(((x-y).array().pow(2)).sum()+(std::pow(cc,2)));
    return(res);
}


double PearsonKernel(Eigen::RowVectorXd x, Eigen::RowVectorXd y, Eigen::RowVectorXd parms)
{
    double constant=parms(0);
    double omega=parms(1);
    
    //Inner product
    double xx = 0;
    double yy = 0;
    double xy = 0;
    for (int i=0;i<x.size();i++)
    {
        double u = x(i) * x(i);
        double v = y(i) * y(i);
        double uv = x(i) * y(i);
        xx += u;
        yy += v;
        xy += uv;
    }
    
    double m = constant * std::sqrt(-2.0 * xy + xx + yy);
    return 1.0 / std::pow(1.0 + m * m, omega);
}


double PolynomialKernel(Eigen::RowVectorXd x, Eigen::RowVectorXd y,Eigen::RowVectorXd parms)
{
    double a=parms(0);
    double b=parms(1);
    double num=x.dot(y);
    double res=std::pow(num+b,a);
    return(res);
}


double PowerKernel(Eigen::RowVectorXd x,Eigen::RowVectorXd y,Eigen::RowVectorXd parms)
{
    double d=parms(0);
    double res = (-1.0*std::pow(std::sqrt(((x-y).array().pow(2)).sum()),d));
    return(res);
}


double RationalQuadraticKernel(Eigen::RowVectorXd x, Eigen::RowVectorXd y, Eigen::RowVectorXd parms)
{
    double cc=parms(0);
    double num=((x-y).array().pow(2)).sum();
    double den=((x-y).array().pow(2)).sum()+cc;
    double res = 1.0-(num/den);
    return(res);
}


double SigmoidKernel(Eigen::RowVectorXd x, Eigen::RowVectorXd y, Eigen::RowVectorXd parms)
{
    double alpha=parms(0);
    double constant=parms(1);
    double sum = 0.0;
    for (int i=0;i<x.size();i++)
        sum += x(i) * y(i);
    double value = std::tanh(alpha * sum + constant);
    
    return value;
}


double SphericalKernel(Eigen::RowVectorXd x,Eigen::RowVectorXd y,Eigen::RowVectorXd parms)
{
    double sigma=parms(0);
    double norm = (x-y).norm();
    double res=0;
    if(norm>=sigma){
        res=1.0-1.5*(norm/sigma)+0.5*std::pow(norm/sigma,3);
    }
    return(res);
}


double SplineKernel(Eigen::RowVectorXd x,Eigen::RowVectorXd y,Eigen::RowVectorXd parms)
{
    double prod=1;
    for(int i=0;i<x.size();i++){
        double res=1+(x(i)*y(i)*std::min(x(i),y(i)))-((x(i)+y(i))/2.0)*(std::pow(std::min(x(i),y(i)),2)+(std::pow(std::min(x(i),y(i)),3)/3.0));
        prod=prod*res;
    }
    return(prod);
}


double SquaredSincKernel(Eigen::RowVectorXd x,Eigen::RowVectorXd y,Eigen::RowVectorXd parms)
{
    double norm = 0.0;
    double gamma=parms(0);
    
    for (int i=0;i<x.size();i++)
    {
        double d = x(i) - y(i);
        norm += d * d;
    }
    
    double num = gamma * std::sqrt(norm);
    double den = gamma * gamma * norm;
    
    return std::sin(num) / den;
}


double SymmetricTriangleKernel(Eigen::RowVectorXd x,Eigen::RowVectorXd y,Eigen::RowVectorXd parms)
{
    double gamma=parms(0);
    double norm = 0.0;
    for (int i=0;i<x.size();i++)
    {
        double d = x(i) - y(i);
        norm += d * d;
    }
    
    double z = 1.0 - gamma * std::sqrt(norm);
    
    return (z > 0) ? z : 0;
}


double ThinSplinePlateKernel(arma::vec x,arma::vec y,arma::vec parms)
{
  double sigma=parms(0);
  double r = 0;
  for (int i=0;i<x.n_elem;i++)
  {
    double dxy = x(i) - y(i);
    r += dxy * dxy;
  }

  return r / (sigma * sigma) * std::log(std::sqrt(r) / sigma);
}


double WaveKernel(arma::vec x,arma::vec y,arma::vec parms)
{
  double norm = arma::norm(x-y);
  double sigma=parms(0);

  if (sigma == 0 || norm == 0)
    return 0;

  return (sigma / norm) * std::sin(norm / sigma);
}


double WaveletKernel(arma::vec x,arma::vec y,arma::vec parms)
{
  double a=parms(0);
  arma::vec hx1=((x-y)/a);
  arma::vec res = arma::cos(1.75*hx1)*arma::exp(-1.0*((arma::square(hx1)/2.0)));
  double res2=1;
  for(int i=0;i<res.n_elem;i++){
    res2=res2*res(i);
  }
  return(res2);
}


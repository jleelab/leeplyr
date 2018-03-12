//for sleep
#include <unistd.h>
//Rcpp
#include <Rcpp.h>
using namespace Rcpp;

#include <cpd/nonrigid.hpp>
//#include <Eigen/Dense>

RcppExport SEXP CoherentPointDriftRegistration(SEXP srcX, SEXP srcY, SEXP dstX, SEXP dstY, SEXP betaR, SEXP lambdaR, SEXP gammaR, SEXP sigmaR, SEXP max_iterR){
   BEGIN_RCPP

   Rcpp::RNGScope __rngScope;

    double beta = Rcpp::as<double>(betaR);
    double lambda = Rcpp::as<double>(lambdaR);
    double gamma = Rcpp::as<double>(gammaR);
    double sigma = Rcpp::as<double>(sigmaR);
    int max_iter = Rcpp::as<int>(max_iterR);
            
  	//resized used for registration
    //int width = Rcpp::as<int>(N); // N
    //int height = Rcpp::as<int>(M); // M

    //define meshgrid in y-coordinates could be done more elegant with Eigen but quick and dirty for now
    /*cpd::Matrix y(height*width, 1);
    cpd::Matrix x(y.rows(),1);
    for (int i = 0; i < width; ++i) {
        for (int j = 0; j < height; ++j) {
            y(j+height*i, 0) = j;
            x(j+height*i, 0) = i;
        }
    }

    //column bind to a 2D matrix
    cpd::Matrix grid(y.rows(), x.cols()+y.cols());
	  grid << x, y;*/


	//reading in the correspondence points
	std::vector<int> srX = as<std::vector<int> >(srcX);
	std::vector<int> srY = as<std::vector<int> >(srcY); 
	std::vector<int> dtX = as<std::vector<int> >(dstX); 
 	std::vector<int> dtY = as<std::vector<int> >(dstY); 


	cpd::Matrix iP(srX.size(), 2),  iiP(srX.size(), 2);
	//iP << srX, srY;
	//iiP << dtX, dtY;
	for (unsigned i=0; i < srX.size(); i++) {
	   iP(i, 0) = srX[i];
	   iP(i, 1) = srY[i];
	   iiP(i, 0) = dtX[i];
	   iiP(i, 1) = dtY[i];   
	}



	//initialize the runner
  	cpd::Nonrigid runner;
  	//set the parameters
//    runner.correspondence(true).outliers(gamma).sigma2(sigma).tolerance(1e-5).max_iterations(max_iter);
    runner.correspondence(true).outliers(gamma).sigma2(sigma).tolerance(1e-5).max_iterations(max_iter);
    cpd::NonrigidResult result = runner.run(iP, iiP);

  	//Rcpp::Rcout << "Results:\n" << result.points << std::endl;
  	//Rcpp::Rcout << "Reference:\n" << m_fish << std::endl;
  	Rcpp::Rcout << "Iterations: " << result.iterations << std::endl;
  	//Rcpp::Rcout << "Correspondance: " << result.correspondence << std::endl;

  	cpd::GaussTransformDirect direct;
    cpd::Probabilities probabilities = direct.compute(iP, result.points, result.sigma2, gamma);

    //Rcpp::Rcout << "Correspondance:\n" << probabilities.correspondence << std::endl;
    //Rcpp::Rcout << "POINTS:\n" << result.points << std::endl;
    

   	Rcpp::Rcout << "\nComputing Coherent Point Shift (CPD) transformation\n" << std::endl;
   	//cpd::Matrix transform = result.transformation_grid(iiP, grid);

   	//from the eigen vector to the std vector
    //std::vector<double> trans(transform.data(), transform.data() + transform.rows() * transform.cols());
    std::vector<int> correspondance(probabilities.correspondence.data(), probabilities.correspondence.data() + probabilities.correspondence.rows() * probabilities.correspondence.cols());
    std::vector<double> points(result.points.data(), result.points.data() + result.points.rows() * result.points.cols());
    
  	return List::create(
  	  _["points"] = points,
    	_["corr.index"] = correspondance
  	);
  	
  	/*   	return List::create(
  	 _["trans"] = trans,
  	_["corr.index"] = correspondance
  	); */

    END_RCPP
    
}
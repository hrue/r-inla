
#include "spatial_covariates.hh"
#include "fmesher.hh"
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <iostream>
#include <fstream>
#include <sstream>
#define WHEREAMI __FILE__ << "(" << __LINE__ << ")\t"
#include <stdio.h>

int main( int argc, char* argv[]) {

  //./covariates mesh_prefix cov_data_prefix

  if (argc > 4){
    std::cout<<"The syntax is ./covariates mesh_prefix conv_data lambda"<<std::endl;
    return 1;
  }
  //  try{
  char* mesh_prefix = argv[1];
  char* fname_prefix = argv[2];
  double lam = 0.0;
  if (argc > 3 ) {
    char* tmp834 = argv[3];
    lam = strtod(tmp834,NULL);
  }


  // }
  /* catch(...){
    std::cout<<"There is a problem here!"<<std::endl;
    return 1;
    }*/
  //read the mesh and make the mesh object
  fmesh::MatrixC matrices;
  matrices.io(true,true);
  matrices.input_prefix(mesh_prefix);
  matrices.load("s");
  matrices.load("tv");
  fmesh::Matrix<double> &S = matrices.DD("s");
  fmesh::Matrix<int>& TV = matrices.DI("tv");



  fmesh::Mesh M(fmesh::Mesh::Mtype_plane,0,true,true);
  M.S_set(S);
  M.TV_set(TV);



  //Read the covariate information.

  fmesh::MatrixC cov_mat;
  cov_mat.io(true,true);
  cov_mat.input_prefix( fname_prefix);
  cov_mat.output_prefix(fname_prefix);
  std::cout<<"fname_prefix = "<<fname_prefix<<std::endl;
  cov_mat.load("dat");
  fmesh::Matrix<double> &tmp = cov_mat.DD("dat");

  std::vector<pointProcess::SpatialCovariate<double> > cov;
  cov.resize(tmp.rows());



  for(int i = 0; i < tmp.rows(); i++) {
    cov[i].set_x(tmp(i,0));
    cov[i].set_y(tmp(i,1));
    cov[i].set_cov(tmp(i,2));
  }

  

  pointProcess::InterpolatedCovariate<double> interp_cov( M,cov);

  //  fmesh::Matrix<double> cov_weights;
  //  interp_cov.weights(cov_weights);

  fmesh::SparseMatrix<double> phiMat;
  interp_cov.phi(phiMat);

  std::cout<<"dim(phiMat) = " <<phiMat.rows()<<", "<<phiMat.cols()<<std::endl;
  fmesh::Matrix<double> missMat;
  interp_cov.outside_points(missMat);

  //std::cout<<WHEREAMI<<"m =  "<< phiMat.rows() <<" n = "<<phiMat.cols() << std::endl;

  // NOT IMPLEMENTED fmesh::SparseMatrix<double> reg;
  //  interp_cov.reg(reg );
  //Output everything.
  std::string cov_name = "cov";

  // fmesh::Matrix<double> temp = cov_weights;
  //  cov_mat.attach(cov_name,&cov_weights,false);
  //  cov_mat.matrixtype(cov_name,fmesh::IOMatrixtype_general);
  //  cov_mat.output(cov_name);
  std::string phi_name = "phi";
  cov_mat.attach(phi_name,&phiMat,false);
  //  cov_mat.matrixtype(cov_name,fmesh::IOMatrixtype_general);
  cov_mat.output(phi_name);

  //std::string reg_name = "Q";
  //cov_mat.attach(reg_name,&reg,false);
  //  cov_mat.matrixtype(cov_name,fmesh::IOMatrixtype_general);
  //cov_mat.output(reg_name);

  std::string miss_name = "miss";
  cov_mat.attach(miss_name, &missMat,false);
  cov_mat.output(miss_name);

  cov_mat.save();

  return 0;
  



    /* matrices.matrixtype(Gname,fmesh::IOMatrixtype_symmetric);
      matrices.output(Gname);*/

  /*//MAKE THIS WORK
 std::vector<pointProcess::SpatialCovariate<double> > cov;
  std::ifstream input(fname);
  std::string lineData;

  while(getline(input,lineData)) {

    pointProcess::SpatialCovariate<double> tmp;
    std::stringstream line(lineData);
    line >> tmp;
    cov.push_back(tmp);
    }*/
}

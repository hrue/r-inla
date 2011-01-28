#ifndef _DAN_COVARIATES_TCC_
#define _DAN_COVARIATES_TCC_ 1
#include "spatial_covariates.hh"
#define WHEREAMI __FILE__ << "(" << __LINE__ << ")\t"
#define _LOG(msg) std::cout << WHEREAMI << msg<<std::endl;

//_LOG("This is what I'm trying to do here." << std::endl)

#define NOT_IMPLEMENTED (std::cout                                      \
                        << __FILE__ << "(" << __LINE__ << ")\t"        \
                        << "NOT IMPLEMENTED: "                         \
                        << __PRETTY_FUNCTION__ << std::endl);


namespace pointProcess {

  template <typename covariateType> InterpolatedCovariate<covariateType>::InterpolatedCovariate ( const fmesh::Mesh& mesh, const std::vector<SpatialCovariate<covariateType> >& covariates) {
    //smooth can be 'W' for white noise, 'D' for de Wijj, 'T' for thin plate spline.
    //Size checks!  l8r
    
    fmesh::Mesh mesh_(mesh);
    std::size_t n_vertex = mesh_.nV(); //number of vertices

    std::cout<<"Debug Text.  spatial_covariates.tcc Line 15. n_vertex = "<<n_vertex<<"."<<std::endl;

    std::size_t size_least_squares = 3*covariates.size(); //each data point will only have at most 3 non-zero entries in least squares matrix.

    std::cout<<"Debug Text.  spatial_covariates.tcc Line 19. size_least_squares = "<<size_least_squares<<"."<<std::endl;

    //Allocate matrices for ij-val storage for least squares problem (convert to crs laer_
    fmesh::Dart location;
    fmesh::Dart start;
    location = fmesh::Dart(mesh_,0);

    //Current hack!

    typedef fmesh::Matrix1<double> vecd;
    typedef fmesh::Matrix1<int> veci;

    veci ivec(size_least_squares);
    veci jvec(size_least_squares);
    vecd valvec(size_least_squares);
    int index = 0;


    std::size_t sights_it=0; //number of empty points outside mesh.
 for(typename std::vector<SpatialCovariate<covariateType> >::const_iterator cov=covariates.begin();cov != covariates.end();cov++) {
     
     


      //Step 1: Find the location of point i

      typename std::vector<SpatialCovariate<covariateType> >::const_iterator::difference_type data_num= distance(covariates.begin(),cov);
      //  std::cout<<"num "<<data_num << ":  ";
      //      covariateVector[data_num] = (*cov).cov_value();
      
      fmesh::Point loc(cov->x(),cov->y(),0.0);
      
      //Find the point using a robust search
      std::size_t search_tri = 0;
      location =fmesh::Dart();
      while( (location.isnull()) && ( search_tri < mesh_.nT()) ) {
	start = fmesh::Dart(mesh_,search_tri);
        location = mesh_.locate_point(start, loc); 
	search_tri++;
      }
      if (!location.isnull()) {
	location = fmesh::Dart(mesh_,location.t());  //check the orientation/ordering

	//Step 2: Get Barycentric coords of the point
	fmesh::Point bary;
	mesh_.barycentric(location,loc,bary);

	//Step 3: Get matrix!
	
	// Get a vertex index: location.t() -= which triangle
	fmesh::Int3 verts = mesh_.TV(location.t());

	for(std::size_t vert_num =0; vert_num <3; vert_num++) {
	  // std::cout<< std::distance(mat_ijval.begin(),ijval_iter)<<std::endl;
	  //Add things appropriately
	  ivec(index) = data_num;
	  jvec(index) = verts[vert_num];
	  valvec(index) = bary[vert_num];
	  //  std::cout<<verts[vert_num]<<",";
	 index++;
	}
	//	std::cout<<std::endl;
      } // end if!
      else {
	empty_sights_(sights_it,1) = data_num;
	sights_it++;
      } //end else
    } //End For


std::cout<<"ivec length = "<<ivec.rows() <<std::endl;

    phiMatrix_.clear().cols(n_vertex).rows(covariates.size()); //////////(covariates.size(),n_vertex); //HERE HERE
    phiMatrix_.fromlist(ivec,jvec,valvec);

std::cout<<"phiMatrix_ dims "<< phiMatrix_.rows() <<" " << phiMatrix_.cols()<<std::endl;

 } //  InterpolatedCovariate( const fmesh::Mesh& mesh, const std::vector<SpatialCovariate>& covariates)

 template<typename covariateType> std::ostream& operator<<  (std::ostream& output,
							       const InterpolatedCovariate<covariateType>& cov) {
   output << cov.mesh_.nv() << std::endl;
   for ( std::size_t i =0; i < cov.rows(); i++){
     output << cov[i] << std::endl;
   }
   return output;

 }


  template<typename covariateType>    bool InterpolatedCovariate<covariateType>::regLSQR_(const typename Eigen::SparseMatrix<covariateType>& A, const typename Eigen::SparseMatrix<covariateType>& Q, const VectorType& b, VectorType& x, const double lambda, const int max_iter, int& n_iter, typename Eigen::NumTraits<covariateType>::Real& residual, const double atol, const double btol, const double conlim) {
    
    if(A.cols() != Q.rows()) {
      std::cout<<"ERROR!!!"<<std::endl;
    }
    
     

    //Do some allocation 
    std::size_t n_reg = Q.rows();
    std::size_t n_mod =b.rows();
    //    std::size_t n = n_reg + n_mod;
    VectorType u_old_mod(n_mod);
    VectorType u_old_reg(n_reg); //make sure this initialises to zero!!!!!!
    VectorType v_old(n_mod);
    VectorType u_mod(n_mod);
    VectorType u_reg(n_reg);
    VectorType v(n_mod);
    VectorType w_old(n_mod);
    VectorType w(n_mod);
    VectorType x_old(n_mod);

    Real beta_old=Real();
    Real beta =Real();
    Real alpha_old=Real();
    Real alpha =Real();
    Real rho =Real();
    //Real rhoBar_old=Real();
    Real rhoBar =Real();
    Real theta =Real();
    Real phi =Real();
    //Real phiBar_old=Real();
    Real phiBar =Real();
    Real s = Real();
    Real c = Real();
    Real norm_A_square = Real();
    Real norm_b = std::sqrt(b.dot(b));
    Real sq_lam = std::sqrt(lambda);

    _LOG("am i here?")

    //Initialise 
    beta_old =norm_b;  //check what the sqrt function is!!
    u_old_mod = b/beta_old;
    //  std::cout<<A.rows()<<" " <<b.rows() <<std::endl;
    v_old.noalias() = A.transpose()*b; // + Q.transpose()*0!
    //std::cout<<0;
    alpha_old = std::sqrt(v_old.dot(v_old));
    v_old /= alpha_old;
    w_old = v_old;

for (int iii = 17;iii<22;iii++) {
	std::cout<<w_old[iii] << " ";
      }
 std::cout<<std::endl;
    x_old = x;  //assuming the copy constructor works!!
    phiBar = beta_old;
    rhoBar = alpha_old;

    n_iter = 0;
    bool terminate = false;
    bool terminate1 = false;
    bool terminate2 = false;
    bool terminate3 = false;


    while ( (n_iter < max_iter) && !terminate) {
      
      //      std::cout<<"iter "<<n_iter<<"  "<<residual<<std::endl;
      //for (int iii = 17;iii<22;iii++) {
      //	std::cout<<x[iii] << "  ";
      //      }
      // std::cout<<std::endl;
      n_iter++;

      //Extend subpsace
      u_mod.noalias() = A*v_old - alpha_old*u_old_mod;
      // std::cout<< 1;
      u_reg.noalias() = sq_lam* (Q*v_old) - alpha_old*u_old_reg;

      // std::cout<< 2;
      beta = std::sqrt(u_mod.dot(u_mod) + u_reg.dot(u_reg)); 

      //  std::cout<<"beta = "<<beta;
      if (beta < 1e-15) {
	return false;
      }
      u_mod = u_mod/beta;
      u_reg = u_reg/beta;
      v.noalias() = A.transpose()*u_mod + sq_lam*(Q.transpose()*u_reg)-beta*v_old;
      alpha = std::sqrt(v.dot(v));
      //std::cout<<", alpha=" << alpha <<std::endl;
      if (alpha <1e-15) { 
	return false;
      }
      v = v/alpha;
      norm_A_square = norm_A_square + alpha_old*alpha_old  + beta*beta;

      //Construct orthogonal transoformation
      rho = std::sqrt(rhoBar*rhoBar + beta*beta);
      c = rhoBar/rho;
      s = beta/rho;
      theta = s*alpha;
      rhoBar = -c*alpha;
      phi = c*phiBar;
      phiBar = s*phiBar;

      //Update x and w;
      x.noalias() = x_old + (phi/rho)*w_old;
      w.noalias() = v - (theta/rho)*w_old;

      //Test for convergence
      terminate1 = (phiBar <= (btol*norm_b + atol*std::sqrt(norm_A_square)*std::sqrt(x.dot(x)))); //make norm x better
       terminate2 = (phiBar*alpha*std::abs(c) <= atol*std::sqrt(norm_A_square)*phiBar);
       terminate3 = false;  //not implemented

      terminate = (terminate1 || terminate2) || terminate3; //  check this carefully!!
      terminate = terminate & (n_iter > 50);

      residual = phiBar;

      //make the new old.
      if (!terminate && (n_iter < max_iter)) {
	x_old=x;
	w_old=w;
	v_old=v;
	u_old_mod=u_mod;
	u_old_reg = u_reg;
      }


    }//End while

    return true;

  }


 template<typename scalar>  static bool compareijval(const ijval<scalar>& w1, const ijval<scalar>& w2) {
  
    if (w1.j < w2.j) {
      return (w1.i < w2.i);
    }

    return false;
  }



  template<typename covariateType> void InterpolatedCovariate<covariateType>::ijvalToEigenSparseMatrix_( typename std::vector<ijval<covariateType> >& mat,  typename Eigen::SparseMatrix<covariateType>& eigenMatrix) {
      eigenMatrix.reserve(mat.size()); //Make sure there's enough memory
      

 std::sort(mat.begin(),mat.end(),compareijval<covariateType>);  //sort by j [for efficiency]

      std::size_t it = 0;


      for (std::size_t j =0; (int) j<eigenMatrix.cols(); j++) {

	eigenMatrix.startVec(j);

	while(mat[it].j == j && it < mat.size()) {
	  std::cout<<"j= "<< j <<" i="<<mat[it].i << std::endl;
	  eigenMatrix.insertBack(mat[it].i,j) = mat[it].val;
	  it++;

	} //end while
      }//end for j
      eigenMatrix.finalize();

    }


  template<typename covariateType> void InterpolatedCovariate<covariateType>::fMesherSparseToEigenSparseMatrix_(const fmesh::SparseMatrix<double>* Qpoint,  typename Eigen::SparseMatrix<covariateType>& QMatrix) {
    fmesh::Matrix1<int> i_fm;
    fmesh::Matrix1<int> j_fm;
    fmesh::Matrix1<double> val_fm;
    int rubbish =  Qpoint->tolist(i_fm,j_fm,val_fm);
    rubbish++;
    typename std::vector<ijval<double> > fm_ijval(i_fm.rows());
    for (int it=0;it < i_fm.rows();it++) {
      fm_ijval[it].i = i_fm[it];
      fm_ijval[it].j  = j_fm[it];
      fm_ijval[it].val = val_fm[it];
    }
    
    QMatrix.reserve(i_fm.rows());

    int it = 0;
    for (std::size_t j =0; (int) j< QMatrix.cols(); j++) {
      QMatrix.startVec(j);
      while((fm_ijval[it].j == j) && (it< i_fm.rows())) {
        QMatrix.insertBack(fm_ijval[it].i,j) = fm_ijval[it].val;
	it++;
      } //end while
    }//end for j
    QMatrix.finalize();

    }


  template<typename covariateType> void InterpolatedCovariate<covariateType>::smoothCovariate( Real lambda,const char smooth) {
    if (true) {
      std::cout<< "This should never happen - smoothCovariate has not been implemented.  Sorry!"<<std::endl;
    }
   //  //Step 5: Get penalty matrix
//     fmesh::SparseMatrix<double> C0, C1, G1, B1;
//     fmesh::Matrix<double> Tareas;
//     mesh_.calcQblocks(C0,C1,G1,B1,Tareas); //THIS IS WASTEFUL
//     std::cout<<"Line 79"<<std::endl;
//     fmesh::SparseMatrix<double>* Qpoint = NULL;

//     fmesh::SparseMatrix<double> tmp; //this will be needed later
//     fmesh::SparseMatrix<double> C0inv ;
//     std::cout<<"Line 90"<<std::endl;
//     C0inv = inverse(C0,true);
//     std::cout<<"Line 92"<<std::endl;
//     tmp = G1*C0inv;
//     std::cout<<"Line 93"<<std::endl;
//     fmesh::SparseMatrix<double>tmptmp =  tmp*G1;
//     switch (smooth) {
//     case 'D':
//       Qpoint = &G1;
//       break;
//     case 'W':
//       Qpoint = &C0;
//       break;
//     case 'T':
      
    
//       Qpoint = &tmptmp;
//       std::cout<<"Line 94"<<std::endl;

      
//       break;
//     default:
//       std::cout<< "You need to enter 'D', 'W' or 'T'.  Defaulting to de Wijj prior ('D')\n";
//       Qpoint = &G1;
//     }//end case
//     _LOG("pont1")
//       Qmat_ = *Qpoint;
//     _LOG("number2")

//       /////
//       ////
//       ////  OLD CODE!! MAKE THIS WORK.
//       ////
//       ////
//       //Keep everything in a useful container that will sort by j
//       typename std::vector<ijval<covariateType> > mat_ijval(size_least_squares);
//     typename std::vector<ijval<covariateType> >::iterator ijval_iter = mat_ijval.begin();

//     //Also put the vectors somewhere useful.
//     typename Eigen::VectorXd covariateVector(covariates.size()); ///BAD

//     std::cout<<"size LS = "<< size_least_squares << " mesh.nV() = " << mesh.nV() <<" covariates.size() = "<<covariates.size() <<std::endl;

//     for(typename std::vector<SpatialCovariate<covariateType> >::const_iterator cov=covariates.begin();cov != covariates.end();cov++) {
     
     


//       //Step 1: Find the location of point i

//       typename std::vector<SpatialCovariate<covariateType> >::const_iterator::difference_type data_num= distance(covariates.begin(),cov);
//       // std::cout<<"num "<<data_num << std::endl;
//       covariateVector[data_num] = (*cov).cov_value();
      
//       fmesh::Point loc(cov->x(),cov->y(),0.0);
      
//       //If the data isn't in random order, probably best to start at last location (maybe??)
//       location = mesh_.locate_point(location, loc); 
//       if (!location.isnull()) {
// 	location = fmesh::Dart(mesh_,location.t());  //check the orientation/ordering

// 	//Step 2: Get Barycentric coords of the point
// 	fmesh::Point bary;
// 	mesh_.barycentric(location,loc,bary);

// 	//Step 3: Get matrix!
	
// 	// Get a vertex index: location.t() -= which triangle
// 	fmesh::Int3 verts = mesh_.TV(location.t());

// 	for(std::size_t vert_num =0; vert_num <3; vert_num++) {
// 	  // std::cout<< std::distance(mat_ijval.begin(),ijval_iter)<<std::endl;
// 	  //Add things appropriately
// 	  ijval_iter->i = data_num;
// 	  ijval_iter->j = verts[vert_num];
// 	  ijval_iter->val = bary[vert_num];
// 	  ijval_iter++;
// 	}
//       } // end if!
//     } //End For
    
//     //Step 4: Re-format (i,j,val) into something more useful!  This should be an efficient method for constructign an Eigen::SparseMatrix.
//     std::cout<<"Line 69"<<std::endl;
    
//     //       typename Eigen::SparseMatrix<covariateType> phiMatrix(covariates.size(),n_vertex);
//     //        std::cout<<"Line 72"<<std::endl;
//     //
//     //	
//     //    ijvalToEigenSparseMatrix_(mat_ijval,phiMatrix);
//     //    std::cout<<"phiMatrix.rows() = "<<phiMatrix.rows()<<", phiMatrix.cols() = "<<phiMatrix.cols()<<std::endl;
//     //   




    
//     DynamicMatrix phiMatrix(covariates.size(),n_vertex);
//     phiMatrix.reserve(size_least_squares);
//     for (typename std::vector<ijval<covariateType> >::iterator it = mat_ijval.begin(); it!=mat_ijval.end();it++) {
//       phiMatrix.coeffRef(it->i,it->j) += it->val;
//     }
//     //Step 5: Get penalty matrix
//     /*   fmesh::SparseMatrix<double> C0, C1, G1, B1;
// 	 fmesh::Matrix<double> Tareas;
// 	 mesh_.calcQblocks(C0,C1,G1,B1,Tareas); //THIS IS WASTEFUL
// 	 std::cout<<"Line 79"<<std::endl;
// 	 fmesh::SparseMatrix<double>* Qpoint = NULL;

// 	 fmesh::SparseMatrix<double> tmp; //this will be needed later
// 	 fmesh::SparseMatrix<double> C0inv ;
// 	 std::cout<<"Line 90"<<std::endl;
// 	 C0inv = inverse(C0,true);
// 	 std::cout<<"Line 92"<<std::endl;
// 	 tmp = G1*C0inv;
// 	 std::cout<<"Line 93"<<std::endl;
// 	 fmesh::SparseMatrix<double>tmptmp =  tmp*G1;
// 	 switch (smooth) {
// 	 case 'D':
// 	 Qpoint = &G1;
// 	 break;
// 	 case 'W':
// 	 Qpoint = &C0;
// 	 break;
// 	 case 'T':
      
    
// 	 Qpoint = &tmptmp;
// 	 std::cout<<"Line 94"<<std::endl;

      
// 	 break;
// 	 default:
// 	 std::cout<< "You need to enter 'D', 'W' or 'T'.  Defaulting to de Wijj prior ('D')\n";
// 	 Qpoint = &G1;
// 	 }//end case
 
// 	 //std::cout<<(*Qpoint) << std::endl;
// 	 */
//     fmesh::Matrix1<int> i_fm;
//     fmesh::Matrix1<int> j_fm;
//     fmesh::Matrix1<double> val_fm;
//     Qpoint->tolist(i_fm,j_fm,val_fm);

//     DynamicMatrix QMatrix(n_vertex,n_vertex);
//     for (int ittt = 0; ittt <i_fm.rows(); ittt++) {
//       QMatrix.coeffRef(i_fm[ittt],j_fm[ittt]) += val_fm[ittt];
//     }

//     // std::cout<<QMatrix<<std::endl;

//     //Solve least squares problem!!
//     //Get a bad estimate of lambda
//     //   lambda = Real();
    
//     //   for (typename std::vector<ijval<covariateType> >::const_iterator it = mat_ijval.begin(); it != mat_ijval.end(); it++) {
//     //    lambda = lambda + (it->val) * (it->val);
//     //  }
//     //   Real tmp = Real();
   

//     // for (int it=0; it < i_fm.rows();it++){
//     //     tmp = tmp + (val_fm[it]) * (val_fm[it]);
//     //   }
//     //    lambda = lambda/tmp;

//     //    std::cout<<lambda<<std::endl;


//     covariateWeights_.resize(n_vertex); //Need to do more than 
//     for (std::size_t i = 0;i<n_vertex;i++){
//       covariateWeights_[i]=0.0;
//     }
//     Real res = -1.0;
//     int n_iter = 0;

//     // lambda=10.0;
//     _LOG("or did i only get here?")
//       //std::cout<<phiMatrix<<std::endl;
//       //covariateWeights_.fill(covariateType());
//       covariateWeights_ = phiMatrix.transpose()*covariateVector;
//     bool success=   regLSQR_(phiMatrix, QMatrix, covariateVector,covariateWeights_, lambda, 5000,n_iter,res, 1e-5, 1e-5);
//     if (!success)
//       {
// 	std::cout << "It didn't work!"<<std::endl;
//       }

   

}


}; /*End namespace pointProcess*/
#endif

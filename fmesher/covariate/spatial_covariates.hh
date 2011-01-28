#ifndef _DAN_COVARIATES_
#define _DAN_COVARIATES_ 1


#define WHEREAMI __FILE__ << "(" << __LINE__ << ")\t"
#include "mesh.hh"
#include "vector.hh"
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <iostream>
#include <iomanip>
#include <vector>

/* Note to myself about R.
R program structure: 
-Get data and covariates. 
-Send off data to make mesh. -- Haavard - is y|x = Ax + b implemented yet?
-Send covariates of to program!
-Get everything back and run!


 */

namespace pointProcess {

  //Forward declarations

  template <typename covariateType> class InterpolatedCovariate;
  template<typename covariateType> std::ostream& operator<<  (std::ostream& output,
				 const InterpolatedCovariate<covariateType>& cov);



  // template<typename covariateType> std::ostream& operator<< <> (std::ostream& output,
  //				 const InterpolatedCovariate<covariateType>& cov);


  template<typename scalar> struct ijval;
  //  template<typename scalar> bool operator<(ijval<scalar> a,ijval<scalar> b);
  template<typename scalar> struct ijval {
    std::size_t i;
    std::size_t j;
    scalar val;
      bool operator< (const ijval<scalar>& other) const  {
	//	std::cout<<"compare <"<<std::endl;
      return j < other.j;
    }
    // inline friend bool (::operator< <>)(ijval<scalar> a,ijval<scalar> b);

  };//ijvalMatrix

 

  


  template<typename covariateType> class SpatialCovariate {
    /*** Class name: SpatialCovariate.
	 Provies a container for covariate data.  The data cannot be modified.
	 
Constructor: template<typename covariateType> SpatialCovariance(double x, double y, covariateType cov);
This takes the two-dimensional coordinates (x,y) and the value of the covariate at that point.

Memeber functions: 
  -double x(); - Returns the x-cordinate
  -double y(); - Returns the y-coordinate
  -covariateType cov(); - Returns the covariate value.

    ***/

    //  friend std::istream& operator>> <>(istream& is, SpatialCovariate<covariateType>& val);
  private:
    
    double x_;
    double y_;
    covariateType cov_;
  public:
    SpatialCovariate():x_(0.0),y_(0.0),cov_(){}
    SpatialCovariate( double x, double y, covariateType cov): x_(x), y_(y) {
      cov_(cov);
    }
  
    void set_x(double x) {x_ =x;}
    void set_y(double y) {y_ = y;}
    void set_cov(covariateType cov){cov_=cov;}

    double x() const { return x_;} 

    double y() const {return y_;}
    double cov_value() const {return cov_;}
     
  };

  /*  template<typename covariateType> std::istream& (::operator>> <>)(istream& is, SpatialCovariate<covariateType>& val){

      };*/

 
  template <typename covariateType>  class InterpolatedCovariate {
  public:
    friend std::ostream& operator<< <> (std::ostream& output,
				 const InterpolatedCovariate<covariateType>& cov);
 //Some convenient typedefs
    typedef  typename Eigen::Matrix<covariateType,Eigen::Dynamic,1> VectorType;//WARNING WARNING!!!
    typedef typename Eigen::NumTraits<covariateType>::Real Real;
    typedef typename Eigen::DynamicSparseMatrix<covariateType> DynamicMatrix;

  private:
   


   //The data
    VectorType covariateWeights_; //covariates[i] ith vertex, 
    fmesh::Mesh mesh_;
    fmesh::SparseMatrix<double> phiMatrix_;
    fmesh::SparseMatrix<double> Qmat_;
    fmesh::Matrix<double> empty_sights_;


    //some private methods: they are convenient.
    void ijvalToEigenSparseMatrix_( typename std::vector<ijval<covariateType> >& mat,  typename Eigen::SparseMatrix<covariateType>& eigenMatrix);
    void fMesherSparseToEigenSparseMatrix_(const fmesh::SparseMatrix<double>* Qpoint,  typename Eigen::SparseMatrix<covariateType>& QMatrix);
    bool regLSQR_(const typename Eigen::SparseMatrix<covariateType>& A, const typename Eigen::SparseMatrix<covariateType>& Q, const VectorType& b, VectorType& x, const double lambda, const int max_iter, int& n_iter, Real& residual, const double atol, const double btol, const double conlim=0.0);

  public:

    //The public methods
    InterpolatedCovariate( const fmesh::Mesh& mesh, const std::vector<SpatialCovariate<covariateType> >& covariates);

    void smoothCovariate( Real lambda, const char smooth);

    void weights(fmesh::Matrix<double>& weight) {
      weight.cols(1);
      weight.rows(covariateWeights_.rows());
      for (int i =0; i < covariateWeights_.rows(); i++){
	weight(i,0) = covariateWeights_[i];
    }

    }

    void phi(fmesh::SparseMatrix<double>& ph) {

      ph = phiMatrix_;
    }

 void reg(fmesh::SparseMatrix<double>& reg) {
   if (true) {
     std::cout<< "This isn't implemented" << std::endl;
     throw;
   }
   reg = Qmat_;
   std::cout<<"here 4 " <<std::endl;
    }

    void outside_points(fmesh::Matrix<double>& miss) {
      miss = empty_sights_;
    }
  };
   


 //  template <typename dataType> class dataMap {
//   private: 
//     fmesh::SparseMatrix<dataType> map_;
//     fmesh::Mesh mesh_;
//   public:
//     dataMap( const fmesh::Mesh& mesh, const std::vector<SpatialCovariate<dataType> > & dat):mesh_(mesh) {
//       int n_data = dat.rows();
//       int n_verts = mesh.nV();
//       map_.clear().cols(n_verts).rows(n_data);

//       fmesh::Matrix1<int> ivec;
//       fmesh::Matrix1<int> jvec;
//       fmesh::Matrix1<dataType> valvec;

//       for(typename std::vector<SpatialCovariate<dataType> >::const_iterator data_it=dat.begin();cov != dat.end();data_it++) {
     
     


//       //Step 1: Find the location of point i

//       typename std::vector<SpatialCovariate<dataType> >::const_iterator::difference_type data_num= distance(dat.begin(),data_it);
//       //  std::cout<<"num "<<data_num << ":  ";
//       //      covariateVector[data_num] = (*cov).cov_value();
      
//       fmesh::Point loc(dat->x(),dat->y(),0.0);
      
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
// 	  ivec(index) = data_num;
// 	  jvec(index) = verts[vert_num];
// 	  valvec(index) = bary[vert_num];
// 	  //  std::cout<<verts[vert_num]<<",";
// 	 index++;
// 	}
// 	//	std::cout<<std::endl;
//       } // end if!
//     } //End For
//  std::cout<<WHEREAMI<<"m =  "<< covariates.size() <<" n = "<<n_vertex << std::endl;

//   map_.fromlist(ivec,jvec,valvec);


  
//     }

//     void getMap(fmesh::SparseMatrix<dataType>& map) {
//       map = _map;
//     }


//   };

  template <typename Real> class quadMatrix {

  private:
    fmesh::Mesh mesh_;
    fmesh::SparseMatrix<Real> quad2_;
   
   
  public:
    quadMatrix(fmesh::Mesh mesh):mesh_(mesh){}
    void quad2(fmesh::SparseMatrix<Real> Q2) { 

      int nT = mesh_.nT();
      Q2 .clear().cols(mesh_.nV()).rows(mesh_.nv());
      int point = 0;
      for (int tri = 0; tri < nT;tri++) {
	fmesh::Int3 row_map = mesh_.TV(tri);
	for (int pnt = 0; pnt < 3; pnt++) {
	  
	}

      }
    }


  };

} /* namespace pointProcess*/ 



#include "spatial_covariates.tcc"
#endif

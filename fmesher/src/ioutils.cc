
#include <cstddef>
#include <cstring>
#include <set>
#include <map>
#include <sstream>
#include <cmath>

#include "vector.hh"
#include "ioutils.hh"

#define WHEREAMI __FILE__ << "(" << __LINE__ << ")\t"

#ifdef DEBUG
#define VECTOR_LOG(msg) std::cout << WHEREAMI << msg;
#else
#define VECTOR_LOG(msg)
#endif


using std::ios;
using std::cout;
using std::cin;
using std::endl;

namespace fmesh {

   // enum Datatype {Datatype_dense=0,
   //                 Datatype_sparse=1};
   //  /*! int/double */
   //  enum Valuetype {Valuetype_int=0,
   //                  Valuetype_double=1};
   //  /*! general/symmentric/diagonal */
   //  enum Matrixtype {Matrixtype_general=0,
   //                   Matrixtype_symmetric=1,
   //                   Matrixtype_diagonal=2};
   //  /*! rowmajor/colmajor */
   //  enum Storagetype {Storagetype_rowmajor=0,
   //                    Storagetype_colmajor=1};


  IOHeader::IOHeader() { def(); }

  IOHeader& IOHeader::def(const int& ref) {
    def();
    valuetype = IOValuetype_int;
    return *this;
  }

  IOHeader& IOHeader::def(const double& ref) {
    def();
    valuetype = IOValuetype_double;
    return *this;
  }


  IOHeader& IOHeader::def() {
    version = IOHEADER_VERSION;
    elems = 0;
    rows = 0;
    cols = 0;
    datatype = -1;
    valuetype = -1;
    matrixtype = -1;
    storagetype = IOStoragetype_rowmajor;
    return *this;
  }


  IOHeader& IOHeader::collection(const MatrixC& C)
  {
    datatype = IODatatype_collection;
    elems = C.output_size();
    rows = -1;
    cols = -1;
    storagetype = -1;
    return *this;
  }



  std::ostream& operator<<(std::ostream& output, const IOHeader& h)
  {
    int header_length = sizeof(IOHeader)/sizeof(int);
    output << header_length << " ";
    const int* ioheader_p = (const int*)&h;
    for (int i=0; i<header_length; i++) {
      output << ioheader_p[i];
      if (i+1<header_length)
	output << " ";
    }
    return output;
  }


  std::istream& operator>>(std::istream& input, IOHeader& h)
  {
    int header_length = sizeof(IOHeader)/sizeof(int);
    int file_header_length;
    input >> file_header_length;
    if (file_header_length < header_length) {
      h.dense(Matrix<int>(0),IOMatrixtype_general);
      int* ioheader_p = (int*)&h;
      for (int i=0; i<file_header_length; i++)
	input >> ioheader_p[i];
    } else {
      int* ioheader_p = (int*)&h;
      for (int i=0; i<header_length; i++)
	input >> ioheader_p[i];
      if (file_header_length > header_length) {
	int* buf = new int[header_length-file_header_length];
	for (int i=0; i<header_length-file_header_length; i++)
	  input >> buf[i];
	delete[] buf;
      }
    }
    return input;
  }




  IOHelperC& IOHelperC::OL(std::ostream& output)
  {
    const IOHeader& h(IOHelper<int>::h_);
    const bool& bin_(IOHelper<int>::binary_);
    if ((h.elems==0) || (!cM_)) {
      return *this;
    }
    for (MatrixC::outputT::const_iterator outi = cM_->output_.begin();
	 outi != cM_->output_.end();
	 ++outi) {
      if (bin_) {
	int string_size = (*outi).length()+1;
	output.write((char*)&string_size, sizeof(string_size));
	output.write((char*)(*outi).c_str(), sizeof(char)*string_size);
      } else {
	output << *outi << std::endl;
      }
    }
    return *this;
  }

  IOHelperC& IOHelperC::IL(std::istream& input)
  {
    const IOHeader& h(IOHelper<int>::h_);
    const bool& bin_(IOHelper<int>::binary_);
    if ((h.elems==0) || (!cM_)) {
      return *this;
    }
    std::string name;
    for (int i=0; i<h.elems; i++) {
      if (bin_) {
	int string_size;
	input.read((char*)&string_size, sizeof(string_size));
	char* buf = new char[string_size];
	input.read(buf, sizeof(char)*string_size);
	list_.push_back(std::string(buf));
	delete[] buf;
      } else {
	input >> name;
	list_.push_back(name);
      }
    }
    return *this;
  }




  IOHelperC& IOHelperC::OD(std::ostream& output)
  {
    const IOHeader& h(IOHelper<int>::h_);
    const bool& bin_(IOHelper<int>::binary_);
    if ((h.elems==0) || (!cM_)) {
      return *this;
    }
    for (MatrixC::outputT::const_iterator outi = cM_->output_.begin();
	 outi != cM_->output_.end();
	 ++outi) {
      const MCC& mcc = *(cM_->coll_.find(*outi)->second);
      if (mcc.info.datatype == IODatatype_dense)
	if (mcc.info.valuetype == IOValuetype_int) {
	  IOHelperM<int> ioh;
	  ioh.cD(&(mcc.DI())).matrixtype(mcc.info.matrixtype);
	  ioh.binary(bin_).OH(output).OD(output);
	} else {
	  IOHelperM<double> ioh;
	  ioh.cD(&mcc.DD()).matrixtype(mcc.info.matrixtype);
	  ioh.binary(bin_).OH(output).OD(output);
	}
      else
	if (mcc.info.valuetype == IOValuetype_int) {
	  IOHelperSM<int> ioh;
	  ioh.cD(&mcc.SI()).matrixtype(mcc.info.matrixtype);
	  ioh.binary(bin_).OH(output).OD(output);
	} else {
	  IOHelperSM<double> ioh;
	  ioh.cD(&mcc.SD()).matrixtype(mcc.info.matrixtype);
	  ioh.binary(bin_).OH(output).OD(output);
	}
    }
    return *this;
  }
  
  IOHelperC& IOHelperC::ID(std::istream& input)
  {
    const IOHeader& h(IOHelper<int>::h_);
    const bool& bin_(IOHelper<int>::binary_);
    if ((!M_) || (h.elems==0)) {
      return *this;
    }
    for (listT::iterator listi = list_.begin();
	 listi != list_.end();
	 ++listi) {

      IOHelper<int> ioh_;
      ioh_.binary(bin_).IH(input);

      if (ioh_.h_.datatype == IODatatype_dense)
	if (ioh_.h_.valuetype == IOValuetype_int) {
	  IOHelperM<int> ioh;
	  ioh.D(&(M_->DI(*listi)));
	  ioh.binary(bin_).IH(ioh_.h_).ID(input);
	} else {
	  IOHelperM<double> ioh;
	  ioh.D(&(M_->DD(*listi)));
	  ioh.binary(bin_).IH(ioh_.h_).ID(input);
	}
      else
	if (ioh_.h_.valuetype == IOValuetype_int) {
	  IOHelperSM<int> ioh;
	  ioh.D(&(M_->SI(*listi)));
	  ioh.binary(bin_).IH(ioh_.h_).ID(input);
	} else {
	  IOHelperSM<double> ioh;
	  ioh.D(&(M_->SD(*listi)));
	  ioh.binary(bin_).IH(ioh_.h_).ID(input);
	}
    }
    return *this;
  }
  




  template <>
  Matrix<int>& MatrixC::attach(std::string name,
			       Matrix<int>* M,
			       bool transfer_ownership,
			       IOMatrixtype matrixt)
  {
    free(name);
    coll_.insert(collPairT(name,
			   new MCC(IODatatype_dense,
				   IOValuetype_int,
				   matrixt,
				   M,transfer_ownership)));
    activate(name);
    return coll_[name]->DI();
  }
  
  template <>
  Matrix<double>& MatrixC::attach(std::string name,
				  Matrix<double>* M,
				  bool transfer_ownership,
				  IOMatrixtype matrixt)
  {
    free(name);
    coll_.insert(collPairT(name,
			   new MCC(IODatatype_dense,
				   IOValuetype_double,
				   matrixt,
				   M,transfer_ownership)));
    activate(name);
    return coll_[name]->DD();
  }



  template <>
  SparseMatrix<int>& MatrixC::attach(std::string name,
				     SparseMatrix<int>* M,
				     bool transfer_ownership,
				     IOMatrixtype matrixt)
  {
    free(name);
    coll_.insert(collPairT(name,
			   new MCC(IODatatype_sparse,
				   IOValuetype_int,
				   matrixt,
				   M,transfer_ownership)));
    activate(name);
    return coll_[name]->SI();
  }
  
  template <>
  SparseMatrix<double>& MatrixC::attach(std::string name,
					SparseMatrix<double>* M,
					bool transfer_ownership,
					IOMatrixtype matrixt)
  {
    free(name);
    coll_.insert(collPairT(name,
			   new MCC(IODatatype_sparse,
				   IOValuetype_double,
				   matrixt,
				   M,transfer_ownership)));
    activate(name);
    return coll_[name]->SD();
  }

  
  bool MatrixC::activate(std::string name)
  {
    collT::iterator colli;
    if ((colli = coll_.find(name)) == coll_.end()) {
      return false;
    }
    colli->second->info.active = true;
    return true;
  }

  void MatrixC::activate()
  {
    for (collT::iterator colli = coll_.begin();
	 colli != coll_.end();
	 ++colli) {
      colli->second->info.active = true;
    }
  }

  void MatrixC::load_file(std::string filename, bool only_list)
  {
    IOHelperC ioh;
    if (filename=="-") {
      /* Can only read stdin once, so read everything now. */
      ioh.D(this);
      ioh.binary(bin_in_).IH(std::cin).IL(std::cin).ID(std::cin);
    } else {
      std::ifstream I;
      I.open(filename.c_str(),
	     (bin_in_ ? (ios::in | ios::binary) : ios::in));
      if (!I.is_open()) {
	// TODO: Add error handling.
      }
      ioh.D(this);
      ioh.binary(bin_in_).IH(I).IL(I);
      if (!only_list)
	ioh.ID(I);
      I.close();
    }

    /* Populate source_ */
    for (IOHelperC::listT::const_iterator listi = ioh.list_.begin();
	 listi != ioh.list_.end();
	 ++listi) {
      source_[(*listi)] = filename;
    }
  }

  MCCInfo MatrixC::load(std::string name)
  {
    /* Is the matrix already loaded? */
    if (activate(name))
      return info(name);

    sourceT::const_iterator sourcei;
    if ((sourcei = source_.find(name)) != source_.end()) {
      /* The matrix is in a collection file */
      load_file(sourcei->second);
      if (activate(name))
	return info(name);
    }

    /* Do we have a prefix to read from? */
    if (input_prefix_ == "-")
      return info(name);

    /* Try to read from prefix data. */
    
    std::ifstream I;
    I.open((input_prefix_+name).c_str(),
	   (bin_in_ ? (ios::in | ios::binary) : ios::in));
    if (!I.is_open()) {
      return info(name);
    }
    IOHelper<int> ioh_;
    ioh_.binary(bin_in_).IH(I);
    I.close();

    coll_.insert(collPairT(name,
			   new MCC(IODatatype(ioh_.h_.datatype),
				   IOValuetype(ioh_.h_.valuetype),
				   IOMatrixtype(ioh_.h_.matrixtype))));
    activate(name);

    I.open((input_prefix_+name).c_str(),
	   (bin_in_ ? (ios::in | ios::binary) : ios::in));
    if (!I.is_open()) {
      return info(name);
    }
    if (ioh_.h_.datatype == IODatatype_dense)
      if (ioh_.h_.valuetype == IOValuetype_int) {
	IOHelperM<int> ioh;
	ioh.D(&DI(name));
	ioh.binary(bin_in_).IH(I);
	ioh.ID(I);
      } else {
	IOHelperM<double> ioh;
	ioh.D(&DD(name));
	ioh.binary(bin_in_).IH(I).ID(I);
      }
    else
      if (ioh_.h_.valuetype == IOValuetype_int) {
	IOHelperSM<int> ioh;
	ioh.D(&SI(name));
	ioh.binary(bin_in_).IH(I).ID(I);
      } else {
	IOHelperSM<double> ioh;
	ioh.D(&SD(name));
	ioh.binary(bin_in_).IH(I).ID(I);
      }
    I.close();

    return info(name);
  }

  MatrixC& MatrixC::free(std::string name)
  {
    dont_output(name);
    
    collT::iterator colli;
    if ((colli = coll_.find(name)) != coll_.end()) {
      delete colli->second;
      coll_.erase(colli);
    }
    return *this;
  }

  MatrixC& MatrixC::dont_output(std::string name)
  {
    outputT::iterator outi;
    if ((outi = output_.find(name)) != output_.end()) {
      output_.erase(outi);
    }
    return *this;
  }

  MatrixC& MatrixC::output(std::string name)
  {
    if (name=="-") {
      output_all_ = true;
      for (collT::iterator colli = coll_.begin();
	   colli != coll_.end();
	   ++colli) {
	if (colli->second->info.active)
	  output_.insert(colli->first);
      }
    } else if (name=="--") {
      output_all_ = true;
      for (collT::iterator colli = coll_.begin();
	   colli != coll_.end();
	   ++colli) {
	if (activate(colli->first))
	  output_.insert(colli->first);
      }
    } else {
      if (info(name).loaded) {
	activate(name);
	if (output_all_) {
	  output_all_ = false;
	  output_.clear();
	}
	output_.insert(name);
      }
    }
    return *this;
  }




  void MatrixC::io(bool bin_in, bool bin_out)
  {
    bin_in_ = bin_in;
    bin_out_ = bin_out;
  }

  void MatrixC::input_prefix(std::string prefix)
  {
    input_prefix_ = prefix;
  }

  void MatrixC::output_prefix(std::string prefix)
  {
    output_prefix_ = prefix;
  }

  void MatrixC::input_file(std::string filename)
  {
    load_file(filename,true);
  }
  void MatrixC::output_file(std::string filename)
  {
    output_file_ = filename;
  }

  void MatrixC::input_raw(std::string name,
			  std::string specification,
			  std::string filename)
  {
    /* Parse raw ascii matrix data and add to collection. */    

    if (specification=="ddgr") {
      Matrix<double>& M = attach(name,new Matrix<double>());
      if (filename=="-")
	input_raw_M(std::cin,M);
      else {
	std::ifstream I;
	I.open(filename.c_str(), std::ios::in);
	if (!I.is_open()) {
	  // TODO: Add error handling.
	}
	input_raw_M(I,M);
	I.close();
      }
    } else if (specification=="digr") {
      Matrix<int>& M = attach(name,new Matrix<int>());
      if (filename=="-")
	input_raw_M(std::cin,M);
      else {
	std::ifstream I;
	I.open(filename.c_str(), std::ios::in);
	if (!I.is_open()) {
	  // TODO: Add error handling.
	}
	input_raw_M(I,M);
	I.close();
      }
    } else
      NOT_IMPLEMENTED;

  }



  void MatrixC::save()
  {
   /* Write the matrix collection to output */
    if (output_prefix_ != "-") {
      for (outputT::const_iterator outi = output_.begin();
	   outi != output_.end();
	   ++outi) {
	MCC& mcc = *(coll_.find(*outi)->second); 
	if (mcc.info.datatype == IODatatype_dense)
	  if (mcc.info.valuetype == IOValuetype_int)
	    save_M((output_prefix_+(*outi)),mcc.DI(),mcc.info,bin_out_);
	  else
	    save_M((output_prefix_+(*outi)),mcc.DD(),mcc.info,bin_out_);
	else
	  if (mcc.info.valuetype == IOValuetype_int)
	    save_SM((output_prefix_+(*outi)),mcc.SI(),mcc.info,bin_out_);
	  else
	    save_SM((output_prefix_+(*outi)),mcc.SD(),mcc.info,bin_out_);
      }
    }
    if (output_file_ != "") {
      if (output_file_ == "-") {
	IOHelperC ioh;
	ioh.cD(this);
	ioh.binary(bin_out_).OH(std::cout).OL(std::cout).OD(std::cout);
      } else {
	std::ofstream O;
	O.open(output_file_.c_str(),
	       (bin_out_ ? (ios::out | ios::binary) : ios::out));
	if (!O.is_open()) {
	  // TODO: Add error handling.
	}
	IOHelperC ioh;
	ioh.cD(this);
	ioh.binary(bin_out_).OH(O).OL(O).OD(O);
	O.close();
      }
    }
  }







  Matrix<int>& MatrixC::DI(std::string name)
  {
    collT::iterator colli;
    if (((colli = coll_.find(name)) != coll_.end()) &&
	(colli->second->info.datatype == IODatatype_dense) &&
	(colli->second->info.valuetype == IOValuetype_int) &&
	(colli->second->info.active)) {
      return colli->second->DI();
    }
    return attach(name,new Matrix<int>());
  }

  Matrix<double>& MatrixC::DD(std::string name)
  {
    collT::iterator colli;
    if (((colli = coll_.find(name)) != coll_.end()) &&
	(colli->second->info.datatype == IODatatype_dense) &&
	(colli->second->info.valuetype == IOValuetype_double) &&
	(colli->second->info.active)) {
      return colli->second->DD();
    }
    return attach(name,new Matrix<double>());
  }

  SparseMatrix<int>& MatrixC::SI(std::string name)
  {
    collT::iterator colli;
    if (((colli = coll_.find(name)) != coll_.end()) &&
	(colli->second->info.datatype == IODatatype_sparse) &&
	(colli->second->info.valuetype == IOValuetype_int) &&
	(colli->second->info.active)) {
      return colli->second->SI();
    }
    return attach(name,new SparseMatrix<int>());
  }

  SparseMatrix<double>& MatrixC::SD(std::string name)
  {
    collT::iterator colli;
    if (((colli = coll_.find(name)) != coll_.end()) &&
	(colli->second->info.datatype == IODatatype_sparse) &&
	(colli->second->info.valuetype == IOValuetype_double) &&
	(colli->second->info.active)) {
      return colli->second->SD();
    }
    return attach(name,new SparseMatrix<double>());
  }


  
  void MatrixC::matrixtype(std::string name, IOMatrixtype matrixt)
  {
    collT::iterator colli;
    if ((colli = coll_.find(name)) != coll_.end())
      colli->second->info.matrixtype = matrixt;
  }

  MCCInfo MatrixC::info(std::string name) const
  {
    collT::const_iterator colli;
    if ((colli = coll_.find(name)) == coll_.end()) {
      return MCCInfo();
    }
    return colli->second->info;
  }






} /* namespace fmesh */

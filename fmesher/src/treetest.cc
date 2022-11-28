#include "trees.h"
#include "locator.h"

using namespace fmesh;

  int main(void) {
    if (false) {
    SBBTree<double> tree(10);
    SBBTree<double>::iterator i = tree.begin();
    SBBTree<double>::const_iterator ci = tree.begin();
    *i = 1.0;
    //    *ci = 2.0;

    typedef std::pair<double,double> segment_type;
    typedef std::vector<segment_type> segment_vector_type;
    typedef std::vector<segment_vector_type> multi_segment_type;
    typedef SegmentSet<double> Simple_type;
    typedef IntervalTree<double> I_type;
    typedef SegmentTree<double, SegmentSet<double> > S_type;
    typedef SegmentTree<double, I_type > SI_type;
    typedef SegmentTree<double, S_type > SS_type;
    typedef SegmentTree<double, SI_type > SSI_type;
    typedef SegmentTree<double, SS_type > SSS_type;
    typedef multi_segment_type bbox_type;
    
    bbox_type bbox(3);
    bbox[0].resize(4);
    bbox[1].resize(4);
    bbox[2].resize(4);
    bbox[0][0] = segment_type(0,5);
    bbox[0][1] = segment_type(0.5,3.2);
    bbox[0][2] = segment_type(2.3,3.4);
    bbox[0][3] = segment_type(0.5,4);
    bbox[1][0] = segment_type(10,10);
    bbox[1][1] = segment_type(11,11);
    bbox[1][2] = segment_type(12,13);
    bbox[1][3] = segment_type(13,14);
    bbox[2][0] = segment_type(0,4);
    bbox[2][1] = segment_type(0,1);
    bbox[2][2] = segment_type(1,2);
    bbox[2][3] = segment_type(2,3);

    typedef SI_type test_type;
    test_type st(bbox.begin());
    std::cout << st;
    st.add_segment(0,4);
    std::cout << st;
    st.build_tree();

    std::cout << st;

    Point s = Point(0.55,0.0,0.0);
    std::vector<double> sv(1);
    sv[0] = s[0];
    {
      for (typename test_type::search_iterator si = st.search(sv.begin()); !si.is_null(); ++si) {
	std::cout << si.is_null() << " " << *si << std::endl;
      }
    }
    
    {
      OrderedSegmentSet<double> oss(bbox.begin());
      oss.add_segment(0,4);
      std::cout << oss;
      oss.build_tree();
      for (typename OrderedSegmentSet<double>::L_search_iterator si = oss.L_search(sv.begin()); !si.is_null(); ++si) {
	std::cout << si.is_null() << " " << *si << std::endl;
      }
      for (typename OrderedSegmentSet<double>::R_search_iterator si = oss.R_search(sv.begin()); !si.is_null(); ++si) {
	std::cout << si.is_null() << " " << *si << std::endl;
      }
    }

  }

    if (true) {
    /*
    I_type it(bbox.begin());
    it.add_segment(0,4);
    it.build_tree();

    std::cout << it;
    */

    //    std::cout << "***********************************" << std::endl;

    {
      Mesh M;
      M.type(fmesh::Mesh::Mtype_sphere);
      //      M.useX11(true,false,500,500,-1.05,1.05,-1.05,1.05);
      M.make_globe(25, 1.0);
      
      std::cout << M << std::endl;
      
      
      int the_dimensions[] = {0,1};
      std::vector<int> dimensions(the_dimensions,
				  the_dimensions +
				  sizeof(the_dimensions) / sizeof(int) );
      

      Point s = Point(1.0,0.0,0.0);
      int t = -1;

      if (true) {
	TriangleLocator locator1(&M, dimensions, true);
	t = locator1.locate(s);
	std::cout << "Point=" << s << std::endl;
	if (t<0) {
	  std::cout << "Triangle not found." << std::endl;
	} else {
	  std::cout << "Triangle #" << t << " (" << M.TV(t)[0] << "," << M.TV(t)[1] << "," << M.TV(t)[2] << ")" << std::endl;
	}
      }

      if (false) {
      TriangleLocator locator2(&M, dimensions, false);
      t = locator2.locate(s);
      std::cout << "Point=" << s << std::endl;
      if (t<0) {
	std::cout << "Triangle not found." << std::endl;
      } else {
	std::cout << "Triangle #" << t << " (" << M.TV(t)[0] << "," << M.TV(t)[1] << "," << M.TV(t)[2] << ")" << std::endl;
      }
      }

      if (true) {
	TriangleLocator locator1(&M, dimensions, true);
	for (int is=0; is < M.nV(); ++is) {
	  t = locator1.locate(M.S(is));
	  //	  std::cout << "Point=" << M.S(is) << std::endl;
	  if (t<0) {
	    std::cout << "Triangle not found." << std::endl;
	  } else {
	    //	    std::cout << "Triangle #" << t << " (" << M.TV(t)[0] << "," << M.TV(t)[1] << "," << M.TV(t)[2] << ")" << std::endl;
	  }
	}
      }


      /*
      std::cout << locator1;
    std::cout << "***********************************" << std::endl;
      std::cout << locator2;
      */
    }
}

  }

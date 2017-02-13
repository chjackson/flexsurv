
#ifndef MAPPLY_4
#define MAPPLY_4

#include <Rcpp.h>

namespace {

  namespace flexsurv {
    
    template <
      int RTYPE_1, bool NA_1, typename T_1,
      int RTYPE_2, bool NA_2, typename T_2,
      int RTYPE_3, bool NA_3, typename T_3,
      int RTYPE_4, bool NA_4, typename T_4,
      typename Function
      >
    class Mapply_4 : public Rcpp::VectorBase<
      Rcpp::traits::r_sexptype_traits<
	typename ::Rcpp::traits::result_of<Function>::type
	>::rtype ,
      true ,
      Mapply_4<RTYPE_1,NA_1,T_1,RTYPE_2,NA_2,T_2,RTYPE_3,NA_3,T_3,RTYPE_4,NA_4,T_4, Function>
      > {
    public:
      typedef typename ::Rcpp::traits::result_of<Function>::type result_type ;
      
      typedef Rcpp::VectorBase<RTYPE_1,NA_1,T_1> VEC_1 ;
      typedef Rcpp::VectorBase<RTYPE_2,NA_2,T_2> VEC_2 ;
      typedef Rcpp::VectorBase<RTYPE_3,NA_3,T_3> VEC_3 ;
      typedef Rcpp::VectorBase<RTYPE_4,NA_4,T_4> VEC_4 ;
      
      typedef typename Rcpp::traits::Extractor<RTYPE_1,NA_1,T_1>::type EXT_1 ;
      typedef typename Rcpp::traits::Extractor<RTYPE_2,NA_2,T_2>::type EXT_2 ;
      typedef typename Rcpp::traits::Extractor<RTYPE_3,NA_3,T_3>::type EXT_3 ;
      typedef typename Rcpp::traits::Extractor<RTYPE_4,NA_4,T_4>::type EXT_4 ;
      
      Mapply_4(const VEC_1& vec_1_,
	       const VEC_2& vec_2_,
	       const VEC_3& vec_3_,
	       const VEC_4& vec_4_,
	       Function fun_ ) :
	vec_1(vec_1_.get_ref()), vec_2(vec_2_.get_ref()),
	vec_3(vec_3_.get_ref()), vec_4(vec_4_.get_ref()),
	fun(fun_){}
      
      inline result_type operator[]( R_xlen_t i ) const {
	return fun( vec_1[i], vec_2[i], vec_3[i], vec_4[i] );
      }
      
      inline R_xlen_t size() const { return vec_1.size() ; }
      
    private:
      const EXT_1& vec_1 ;
      const EXT_2& vec_2 ;
      const EXT_3& vec_3 ;
      const EXT_4& vec_4 ;
      Function fun ;
    } ;
  }
  
  template <
    int RTYPE_1, bool NA_1, typename T_1,
    int RTYPE_2, bool NA_2, typename T_2,
    int RTYPE_3, bool NA_3, typename T_3,
    int RTYPE_4, bool NA_4, typename T_4,
    typename Function>
  inline
  flexsurv::Mapply_4<RTYPE_1,NA_1,T_1,RTYPE_2,NA_2,T_2,RTYPE_3,NA_3,T_3,RTYPE_4,NA_4,T_4,Function>
  mapply(
	 const Rcpp::VectorBase<RTYPE_1,NA_1,T_1>& t1,
	 const Rcpp::VectorBase<RTYPE_2,NA_2,T_2>& t2,
	 const Rcpp::VectorBase<RTYPE_3,NA_3,T_3>& t3,
	 const Rcpp::VectorBase<RTYPE_4,NA_4,T_4>& t4,
	 Function fun
	 ){
    return flexsurv::Mapply_4<RTYPE_1,NA_1,T_1,RTYPE_2,NA_2,T_2,RTYPE_3,NA_3,T_3,RTYPE_4,NA_4,T_4,Function>( t1, t2, t3, t4, fun ) ;
  }
  
}

#endif

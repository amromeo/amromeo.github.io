// SPDX-License-Identifier: GPL-3.0-only
// This addition is entered only after a repeated Arm A line-search failure.
#define round3_native r3_standard_native
#include "native.cpp"
#undef round3_native
#include <sstream>
#include <iomanip>

namespace {
struct Affine {
    Objective &f;std::vector<double> center,scale,work;
    void eval(const double *t){
        for(size_t j=0;j<work.size();++j)work[j]=center[j]+scale[j]*t[j];
        f.eval(work.data());
    }
};
double afn(int,double *t,void *v){auto &a=*static_cast<Affine*>(v);a.eval(t);return a.f.value;}
void agr(int n,double *t,double *out,void *v){auto &a=*static_cast<Affine*>(v);a.eval(t);
    for(int j=0;j<n;++j)out[j]=a.f.gradient[j]*a.scale[j];}
}
extern "C" SEXP round3_native(SEXP xs,SEXP ms,SEXP as,SEXP spec,SEXP start,SEXP lo,SEXP hi,
                              SEXP rule,SEXP max_iterations,SEXP run_mode) {
    SEXP first=PROTECT(r3_standard_native(xs,ms,as,spec,start,lo,hi,rule,max_iterations,run_mode));
    if(Rf_length(rule)!=0||INTEGER(run_mode)[0]!=1||INTEGER(max_iterations)[0]!=10000||
       INTEGER(VECTOR_ELT(first,3))[0]!=52){UNPROTECT(1);return first;}
    SEXP second=PROTECT(r3_standard_native(xs,ms,as,spec,VECTOR_ELT(first,0),lo,hi,rule,max_iterations,run_mode));
    const double initial=REAL(VECTOR_ELT(first,1))[0],ll1=REAL(VECTOR_ELT(first,2))[0],ll2=REAL(VECTOR_ELT(second,2))[0];
    const int code2=INTEGER(VECTOR_ELT(second,3))[0];
    std::ostringstream trace;trace<<std::setprecision(17)<<"ordinary code 52, logLik "<<ll1<<", retry code "<<code2<<", logLik "<<ll2;
    SET_VECTOR_ELT(second,1,Rf_ScalarReal(initial));
    SET_VECTOR_ELT(second,6,VECTOR_ELT(first,6));SET_VECTOR_ELT(second,8,VECTOR_ELT(first,8));
    if(code2!=52) {
        trace<<", "<<CHAR(STRING_ELT(VECTOR_ELT(second,4),0));
        SET_VECTOR_ELT(second,4,Rf_mkString(trace.str().c_str()));UNPROTECT(2);return second;
    }
    const double *s=REAL(spec);int n=Rf_length(start),k=(int)s[0],g=(int)s[1];
    Objective f{REAL(xs),REAL(ms),REAL(as),Rf_length(xs),k,g,n,(int)s[2],
                s[3],s[4],s[5],s[6],s[7],s[8],bool(ISNAN(s[7])),bool(ISNAN(s[8])),{},{},0};
    const double *p=REAL(VECTOR_ELT(second,0));
    Affine a{f,{p,p+n},std::vector<double>(n,1),std::vector<double>(n)};
    std::vector<double> t(n,0),lower(n),upper(n),left(p,p+n),right(p,p+n);std::vector<int> box(n,2);
    const double step=std::cbrt(std::numeric_limits<double>::epsilon());
    for(int j=0;j<n;++j) {
        const double h=step*std::max(1.0,std::abs(p[j]));
        left[j]=std::max(REAL(lo)[j],p[j]-h);right[j]=std::min(REAL(hi)[j],p[j]+h);
        if(right[j]>left[j]) {
            f.eval(left.data());const double low=f.gradient[j];f.eval(right.data());
            const double curv=(f.gradient[j]-low)/(right[j]-left[j]);
            if(std::isfinite(curv)&&curv>0)a.scale[j]=1/std::sqrt(curv);
        }
        left[j]=right[j]=p[j];lower[j]=(REAL(lo)[j]-p[j])/a.scale[j];upper[j]=(REAL(hi)[j]-p[j])/a.scale[j];
    }
    double value;int code=0,fn=0,gr=0;char msg[160]="";
    lbfgsb(n,5,t.data(),lower.data(),upper.data(),box.data(),&value,afn,agr,&code,&a,1000,1e-9,&fn,&gr,10000,msg,0,10);
    a.eval(t.data());trace<<", affine code "<<code<<", logLik "<<-f.value*f.n<<", "<<msg;
    SET_VECTOR_ELT(second,0,as_vec(a.work));SET_VECTOR_ELT(second,2,Rf_ScalarReal(-f.value*f.n));
    SET_VECTOR_ELT(second,3,Rf_ScalarInteger(code));SET_VECTOR_ELT(second,4,Rf_mkString(trace.str().c_str()));
    SET_VECTOR_ELT(second,5,as_vec(f.gradient));SET_VECTOR_ELT(second,7,as_vec(a.work));
    SET_VECTOR_ELT(second,9,Rf_ScalarInteger(fn+INTEGER(VECTOR_ELT(first,9))[0]+INTEGER(VECTOR_ELT(second,9))[0]));
    UNPROTECT(2);return second;
}

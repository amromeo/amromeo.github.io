// SPDX-License-Identifier: GPL-3.0-only
// Original-scale anchored normal-mixture log likelihood and analytic derivative.
// Compress routine and anchor counts separately. Component one receives every anchor.
#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Applic.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>
#include <cstring>

// The archived Arm B build disabled fused multiply-add contraction.
#if defined(__clang__)
#pragma clang fp contract(off)
#elif defined(__GNUC__)
#pragma GCC optimize ("fp-contract=off")
#endif

namespace {
struct Objective {
    const double *x, *mass, *anchor;
    int nx, k, groups, npar, n;
    double center, scale, location, spread, fixed_lambda, fixed_shift;
    bool shape, shifted;
    std::vector<double> last, gradient;
    double value;
    bool cached = false;

    void eval(const double *par) {
        if (cached && std::equal(last.begin(), last.end(), par)) return;
        last.assign(par, par+npar);
        cached = true;
        gradient.assign(npar, 0.0);
        const int shape_i = 2*k + groups*(k-1);
        const double lam = shape ? par[shape_i] : fixed_lambda;
        const double d = shifted ? location + spread*par[npar-1] : fixed_shift;
        const bool linear = !shape && lam == 1.0;
        const double lnscale = std::log(scale);
        constexpr double log_root_two_pi = 0.91893853320467274178;
        double invsd[5], weights[10], logw[10], group_counts[2] = {};
        double ll = 0;
        for (int j=0;j<k;++j) invsd[j] = std::exp(-par[k+j]);
        for (int h=0;h<groups;++h) {
            double sum = 0;
            for (int j=0;j<k;++j) {
                double eta = j==k-1 ? 0 : par[2*k+h*(k-1)+j];
                logw[h*k+j] = eta;
                sum += std::exp(eta);
            }
            for (int j=0;j<k;++j) {
                logw[h*k+j] -= std::log(sum);
                weights[h*k+j] = std::exp(logw[h*k+j]);
            }
        }
        for (int i=0;i<nx;++i) {
            double t=x[i]-d, lt=0, yraw, dl=0, dt=0;
            if (!linear && !(t>0)) { value=1e100; std::fill(gradient.begin(),gradient.end(),0.0); return; }
            if (linear) yraw=t-1;
            else {
                lt=std::log(t);
                const double v=lam*lt;
                if (std::abs(v)<1e-5) {
                    yraw=lt*(1+v*(.5+v*(1.0/6+v*(1.0/24+v/120))));
                    if(shape) dl=lt*lt*(.5+v*(1.0/3+v*(.125+v/30)));
                } else {
                    const double em=std::expm1(v);
                    yraw=lam==1 ? t-1 : em/lam;
                    if(shape) dl=(v*(em+1)-em)/(lam*lam);
                }
                if(shifted) dt=std::exp((lam-1)*lt);
            }
            double y=(yraw-center)/scale, u[5], ld[5], mix[5]={}, count=anchor[i];
            for (int j=0;j<k;++j) {
                u[j]=(y-par[j])*invsd[j];
                ld[j]=-log_root_two_pi-par[k+j]-.5*u[j]*u[j];
            }
            if(anchor[i]>0) {ll+=anchor[i]*ld[0]; mix[0]=anchor[i];}
            for (int h=0;h<groups;++h) {
                const double obs=mass[i+nx*h];
                if(obs==0) continue;
                group_counts[h]+=obs; count+=obs;
                double v[5], peak=-std::numeric_limits<double>::infinity(), sum=0;
                for(int j=0;j<k;++j) { v[j]=ld[j]+logw[h*k+j]; peak=std::max(peak,v[j]); }
                for(int j=0;j<k;++j) { v[j]=std::exp(v[j]-peak); sum+=v[j]; }
                ll+=obs*(peak+std::log(sum));
                for(int j=0;j<k;++j) {
                    const double resp=obs*v[j]/sum;
                    mix[j]+=resp;
                    if(j<k-1) gradient[2*k+h*(k-1)+j]+=resp;
                }
            }
            ll+=count*((lam-1)*lt-lnscale);
            double yscore=0;
            for(int j=0;j<k;++j) {
                const double score=mix[j]*u[j]*invsd[j];
                gradient[j]+=score; gradient[k+j]+=mix[j]*(u[j]*u[j]-1);
                yscore-=score/scale;
            }
            if(shape) gradient[shape_i]+=yscore*dl+count*lt;
            if(shifted) gradient[npar-1]+=(-yscore*dt-count*(lam-1)/t)*spread;
        }
        for(int h=0;h<groups;++h) for(int j=0;j<k-1;++j)
            gradient[2*k+h*(k-1)+j]-=group_counts[h]*weights[h*k+j];
        value=-ll/n;
        for(double &v:gradient) v=-v/n;
        if(!std::isfinite(value)) {value=1e100;std::fill(gradient.begin(),gradient.end(),0.0);}
    }
};
struct Floor {
    double log_min, mean_derivative, shift_derivative;
};

struct Round3 {
    Objective base;
    const double *lo, *hi, *rule;
    bool constrained;
    std::vector<double> current, physical, grad;
    double value=0;
    bool feasible=true;

    Floor floor_at(const double *p, int j) const {
        const double mu=base.center+base.scale*p[j];
        const bool identity=!base.shape && base.fixed_lambda==1;
        const double shift=base.shifted ? base.location+base.spread*p[base.npar-1] : 0;
        const double u=identity ? 0 : std::exp(mu);
        const double median=identity ? mu+1 : shift+u;
        // rule supplies delta, the old identity SD floor and C(b) in the protocol.
        const double fixed_min=std::max(rule[0],rule[1]);
        const double target=std::max(fixed_min,std::max(0.0,median)*rule[2]);
        const double slope=median*rule[2]>fixed_min ? rule[2]/target : 0;
        Floor f{};
        if(identity) f={std::log(target/base.scale),slope*base.scale,0};
        else {
            const double ratio=target/u;
            const double v=2*ratio*ratio/(1+std::hypot(1.0,2*ratio));
            const double log_y=std::log1p(v);
            const double elastic=v==0 ? 1 : v/((1+2*v)*log_y);
            f={.5*std::log(log_y)-std::log(base.scale),elastic*(slope*u-1)*base.scale,
               base.shifted ? elastic*slope*base.spread : 0};
        }
        if(f.log_min<=lo[base.k+j]) return {lo[base.k+j],0,0};
        return f;
    }
    void eval(const double *q) {
        if(!current.empty() && std::equal(current.begin(),current.end(),q))return;
        current.assign(q,q+base.npar);physical=current;
        Floor floors[5]; double penalty=0;
        grad.assign(base.npar,0);
        if(constrained) for(int j=0;j<base.k;++j) {
            floors[j]=floor_at(q,j);
            const int h=base.k+j;
            const double span=hi[h]-floors[j].log_min;
            if(span<0) {
                penalty+=span*span;
                grad[j]+=-2e4*span*floors[j].mean_derivative;
                if(base.shifted)grad.back()+=-2e4*span*floors[j].shift_derivative;
            }
            physical[h]=floors[j].log_min+q[h]*std::max(0.0,span);
        }
        feasible=penalty==0;
        if(!feasible){value=1e6+1e4*penalty;return;}
        base.eval(physical.data());value=base.value;grad=base.gradient;
        if(constrained)for(int j=0;j<base.k;++j) {
            int h=base.k+j;const double old=grad[h];const Floor f=floors[j];
            grad[h]=old*(hi[h]-f.log_min);
            grad[j]+=old*(1-q[h])*f.mean_derivative;
            if(base.shifted)grad.back()+=old*(1-q[h])*f.shift_derivative;
        }
    }
};
double r3_fn(int,double *p,void *v){Round3 &f=*static_cast<Round3*>(v);f.eval(p);return f.value;}
void r3_gr(int n,double *p,double *out,void *v){Round3 &f=*static_cast<Round3*>(v);f.eval(p);std::copy(f.grad.begin(),f.grad.begin()+n,out);}
SEXP as_vec(const std::vector<double> &x){
    SEXP v=Rf_allocVector(REALSXP,x.size());std::copy(x.begin(),x.end(),REAL(v));return v;
}
}

// Mode 0 evaluates a physical start, 1 optimizes it, and 2 evaluates optimizer coordinates.
extern "C" SEXP anchored_fit_start(SEXP xs,SEXP ms,SEXP as,SEXP specs,SEXP start,SEXP lower,SEXP upper,
                              SEXP rules,SEXP max_iterations,SEXP run_mode) {
    if (TYPEOF(specs) != REALSXP || Rf_length(specs) != 9 ||
        TYPEOF(start) != REALSXP || TYPEOF(lower) != REALSXP || TYPEOF(upper) != REALSXP ||
        TYPEOF(xs) != REALSXP || TYPEOF(ms) != REALSXP || TYPEOF(as) != REALSXP ||
        TYPEOF(rules) != REALSXP || (Rf_length(rules) != 0 && Rf_length(rules) != 3) ||
        TYPEOF(max_iterations) != INTSXP || Rf_length(max_iterations) != 1 ||
        TYPEOF(run_mode) != INTSXP || Rf_length(run_mode) != 1)
        Rf_error("Invalid arguments to the anchored mixture optimizer.");
    const double *s=REAL(specs);int k=(int)s[0],g=(int)s[1],np=Rf_length(start),mode=INTEGER(run_mode)[0];
    if(k<1||k>5||g<1||g>2||Rf_length(ms)!=g*Rf_length(xs)||Rf_length(as)!=Rf_length(xs)||
       np!=2*k+g*(k-1)+ISNAN(s[7])+ISNAN(s[8])||Rf_length(lower)!=np||Rf_length(upper)!=np)
        Rf_error("Invalid anchored mixture dimensions.");
    bool constrained=Rf_length(rules)==3;
    if(constrained&&ISNAN(s[7]))Rf_error("Rounding floors require a fixed transform.");
    Round3 f{{REAL(xs),REAL(ms),REAL(as),Rf_length(xs),k,g,np,(int)s[2],s[3],s[4],s[5],s[6],s[7],s[8],
              bool(ISNAN(s[7])),bool(ISNAN(s[8])),{},{},0},REAL(lower),REAL(upper),
              constrained?REAL(rules):nullptr,constrained,{},{},{}};
    std::vector<double> p(REAL(start),REAL(start)+np),lo(REAL(lower),REAL(lower)+np),hi(REAL(upper),REAL(upper)+np);
    std::vector<int> nbd(np,2);
    if(mode!=2)for(int j=0;j<np;++j)p[j]=std::max(lo[j],std::min(hi[j],p[j]));
    if(constrained)for(int j=0;j<k;++j) {
        int h=k+j;double floor=f.floor_at(p.data(),j).log_min;
        if(mode!=2)p[h]=hi[h]>floor?std::clamp((p[h]-floor)/(hi[h]-floor),0.0,1.0):0;
        lo[h]=0;hi[h]=1;
    }
    const auto initial_q=p;f.eval(p.data());const auto initial_physical=f.physical;
    double initial_ll=f.feasible?-f.value*f.base.n:R_NegInf; double result=f.value;
    int code=0,fn=0,gr=0;char msg[160]="evaluation only";
    if(mode==1) {
        lbfgsb(np,5,p.data(),lo.data(),hi.data(),nbd.data(),&result,r3_fn,r3_gr,&code,&f,
               1000.0,1e-9,&fn,&gr,INTEGER(max_iterations)[0],msg,0,10);
        f.eval(p.data());
    }
    SEXP ret=PROTECT(Rf_allocVector(VECSXP,10));
    SET_VECTOR_ELT(ret,0,as_vec(f.physical));SET_VECTOR_ELT(ret,1,Rf_ScalarReal(initial_ll));
    SET_VECTOR_ELT(ret,2,Rf_ScalarReal(f.feasible?-f.value*f.base.n:R_NegInf));
    SET_VECTOR_ELT(ret,3,Rf_ScalarInteger(code));SET_VECTOR_ELT(ret,4,Rf_mkString(msg));
    SET_VECTOR_ELT(ret,5,as_vec(f.grad));SET_VECTOR_ELT(ret,6,as_vec(initial_physical));
    SET_VECTOR_ELT(ret,7,as_vec(p));SET_VECTOR_ELT(ret,8,as_vec(initial_q));SET_VECTOR_ELT(ret,9,Rf_ScalarInteger(fn));
    UNPROTECT(1);return ret;
}

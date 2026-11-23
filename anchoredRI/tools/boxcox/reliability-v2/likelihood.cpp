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
double objective(int, double *p, void *v) {
    Objective &f=*static_cast<Objective*>(v); f.eval(p); return f.value;
}
void derivative(int n, double *p, double *out, void *v) {
    Objective &f=*static_cast<Objective*>(v); f.eval(p);
    std::copy(f.gradient.begin(),f.gradient.begin()+n,out);
}
}

// spec = K, G, n, center, scale, shift location, shift scale, lambda, shift
extern "C" SEXP reliability_fit(SEXP xs, SEXP ms, SEXP as, SEXP specs, SEXP start,
                                SEXP lower, SEXP upper, SEXP optimize) {
    const double *s=REAL(specs);
    const int np=Rf_length(start),k=static_cast<int>(s[0]),g=static_cast<int>(s[1]);
    if(k<1||k>5||g<1||g>2||Rf_length(ms)!=Rf_length(xs)*g ||
       Rf_length(as)!=Rf_length(xs)||np!=2*k+g*(k-1)+ISNAN(s[7])+ISNAN(s[8]))
        Rf_error("Invalid reliability-v2 dimensions");
    Objective f{REAL(xs),REAL(ms),REAL(as),Rf_length(xs),k,g,np,static_cast<int>(s[2]),
                s[3],s[4],s[5],s[6],s[7],s[8],bool(ISNAN(s[7])),bool(ISNAN(s[8])),{},{},0};
    std::vector<double> p(REAL(start),REAL(start)+np);
    std::vector<int> nbd(np,2);
    f.eval(p.data());
    double initial=f.value, final=initial;
    int code=0,fn=0,gr=0; char msg[128]="evaluation only";
    if(LOGICAL(optimize)[0]) {
        lbfgsb(np,5,p.data(),REAL(lower),REAL(upper),nbd.data(),&final,objective,derivative,
               &code,&f,1000.0,1e-9,&fn,&gr,1000,msg,0,10);
        f.eval(p.data()); final=f.value;
    }
    SEXP ans=PROTECT(Rf_allocVector(VECSXP,6));
    SEXP params=PROTECT(Rf_allocVector(REALSXP,np));
    std::copy(p.begin(),p.end(),REAL(params)); SET_VECTOR_ELT(ans,0,params);
    SET_VECTOR_ELT(ans,1,Rf_ScalarReal(-initial*f.n)); SET_VECTOR_ELT(ans,2,Rf_ScalarReal(-final*f.n));
    SET_VECTOR_ELT(ans,3,Rf_ScalarInteger(code));SET_VECTOR_ELT(ans,4,Rf_mkString(msg));
    SEXP gradient=PROTECT(Rf_allocVector(REALSXP,np));
    std::copy(f.gradient.begin(),f.gradient.end(),REAL(gradient));SET_VECTOR_ELT(ans,5,gradient);
    UNPROTECT(3); return ans;
}

//modified zbrent to match needs (need to include params only)
double zbrentNEW_eq(double func(const double,void *),
             const double x1, const double x2,  double f11, double f22, const int compute_first, const double tol,
             void * params)
{
    const int ITMAX=100;
    const double EPS=0.000000001;
    int iter;
    double a=x1,b=x2,c=x2,d,e,min1,min2;
    double fc,p,q,r,s,tol1,xm;
    
    if(compute_first == 0){
        double fa=func(a,params),fb=func(b,params);
    }
    if(compute_first == 1){
        double fa = f11;
        double fb = func(b,params);
    }
    if(compute_first == 2){
        double fa = func(a,params);
        double fb = f22;
    }
    
    
    //if ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0))
    //    {
    //    printf("configucur %d\t borne min : %f\t borne max: %f\n",confi,x1,x2);
    //    nrerror("Root must be bracketed in zbrent");
    //    }
//    if (fa > 0.0 && fb > 0.0)
//    {
//        printf("fa(min) = %f, fb(max) = %f, fa+0.0001 %f, fb-0.0001 %f, borne min : %f\t borne max: %f\n",fa, fb, func(a+0.0001,params), func(b-0.0001,params), x1,x2);
//        NR::nrerror("Root must be bracketed in zbrent 1");
//    }
//    if (fa < 0.0 && fb < 0.0)
//    {
//        printf("fa(min) = %f, fb(max) = %f, fa+0.0001 %f, fb-0.0001 %f, borne min : %f\t borne max: %f\n",fa, fb, func(a+0.0001,params), func(b-0.0001,params),x1,x2);
//        NR::nrerror("Root must be bracketed in zbrent 2");
//    }
    fc=fb;
    for (iter=0;iter<ITMAX;iter++) {
        if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0)) {
            c=a;
            fc=fa;
            e=d=b-a;
        }
        if (fabs(fc) < fabs(fb)) {
            a=b;
            b=c;
            c=a;
            fa=fb;
            fb=fc;
            fc=fa;
        }
        tol1=2.0*EPS*fabs(b)+0.5*tol;
        xm=0.5*(c-b);
        if (fabs(xm) <= tol1 || fb == 0.0) return b;
        if (fabs(e) >= tol1 && fabs(fa) > fabs(fb)) {
            s=fb/fa;
            if (a == c) {
                p=2.0*xm*s;
                q=1.0-s;
            } else {
                q=fa/fc;
                r=fb/fc;
                p=s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0));
                q=(q-1.0)*(r-1.0)*(s-1.0);
            }
            if (p > 0.0) q = -q;
            p=fabs(p);
            min1=3.0*xm*q-fabs(tol1*q);
            min2=fabs(e*q);
            if (2.0*p < (min1 < min2 ? min1 : min2)) {
                e=d;
                d=p/q;
            } else {
                d=xm;
                e=d;
            }
        } else {
            d=xm;
            e=d;
        }
        a=b;
        fa=fb;
        if (fabs(d) > tol1)
            b += d;
        else
            b += SIGN(tol1,xm);
        fb=func(b,params);
    }
//    NR::nrerror("Maximum number of iterations exceeded in zbrent");
    return 0.0;
}





double zbrentNEW(double func(const double,void *),
             const double x1, const double x2, const double tol,
             void * params)
{
    const int ITMAX=100;
    const double EPS=0.000000001;
    int iter;
    double a=x1,b=x2,c=x2,d,e,min1,min2;
    double fa=func(a,params),fb=func(b,params),fc,p,q,r,s,tol1,xm;
    
    //if ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0))
    //    {
    //    printf("configucur %d\t borne min : %f\t borne max: %f\n",confi,x1,x2);
    //    nrerror("Root must be bracketed in zbrent");
    //    }
//    if (fa > 0.0 && fb > 0.0)
//    {
//        printf("fa(min) = %f, fb(max) = %f, fa+0.0001 %f, fb-0.0001 %f, borne min : %f\t borne max: %f\n",fa, fb, func(a+0.0001,params), func(b-0.0001,params), x1,x2);
//        NR::nrerror("Root must be bracketed in zbrent 1");
//    }
//    if (fa < 0.0 && fb < 0.0)
//    {
//        printf("fa(min) = %f, fb(max) = %f, fa+0.0001 %f, fb-0.0001 %f, borne min : %f\t borne max: %f\n",fa, fb, func(a+0.0001,params), func(b-0.0001,params),x1,x2);
//        NR::nrerror("Root must be bracketed in zbrent 2");
//    }
    fc=fb;
    for (iter=0;iter<ITMAX;iter++) {
        if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0)) {
            c=a;
            fc=fa;
            e=d=b-a;
        }
        if (fabs(fc) < fabs(fb)) {
            a=b;
            b=c;
            c=a;
            fa=fb;
            fb=fc;
            fc=fa;
        }
        tol1=2.0*EPS*fabs(b)+0.5*tol;
        xm=0.5*(c-b);
        if (fabs(xm) <= tol1 || fb == 0.0) return b;
        if (fabs(e) >= tol1 && fabs(fa) > fabs(fb)) {
            s=fb/fa;
            if (a == c) {
                p=2.0*xm*s;
                q=1.0-s;
            } else {
                q=fa/fc;
                r=fb/fc;
                p=s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0));
                q=(q-1.0)*(r-1.0)*(s-1.0);
            }
            if (p > 0.0) q = -q;
            p=fabs(p);
            min1=3.0*xm*q-fabs(tol1*q);
            min2=fabs(e*q);
            if (2.0*p < (min1 < min2 ? min1 : min2)) {
                e=d;
                d=p/q;
            } else {
                d=xm;
                e=d;
            }
        } else {
            d=xm;
            e=d;
        }
        a=b;
        fa=fb;
        if (fabs(d) > tol1)
            b += d;
        else
            b += SIGN(tol1,xm);
        fb=func(b,params);
    }
//    NR::nrerror("Maximum number of iterations exceeded in zbrent");
    return 0.0;
}



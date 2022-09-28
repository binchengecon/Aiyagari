
// modified zbrent to match needs (need to include params only)
// Zbrent on N independent conditions.

void zbrentNEW_ND(double func(const double,const double,void *),
             const double x1[l_Z], const double x2[l_Z], double sol[l_Z],
             const double tol,
             void * params)
{

    const int ITMAX=100;
    const double EPS=0.00000000001;
    
    int iter;
    
    double fc,p,q,r,s,tol1,xm,d,e,min1,min2;
    
    double a[l_Z],b[l_Z],c[l_Z];
    double fa[l_Z],fb[l_Z],c[l_Z];
    
    
    func(a,fa,params);
    func(b,fb,params);
    func(c,fc,params);
    
    double e1,d1;
    
    fc1 = fb1;
    fc2 = fb2;
    
    // check that interval is not too narrow?
    
    for (iter=0;iter<ITMAX;iter++) {
    
        for(int Z =0; Z<l_Z; Z++){
            if ((fb[Z] > 0.0 && fc[Z] > 0.0) || (fb[Z] < 0.0 && fc[Z] < 0.0)) {
                c[Z]=a[Z];
                fc[Z]=fa[Z];
                e[Z]=d[Z]=b[Z]-a[Z];
            }
            if (fabs(fc[Z]) < fabs(fb[Z])) {
                a[Z]=b[Z];
                b[Z]=c[Z];
                c[Z]=a[Z];
                fa[Z]=fb[Z];
                fb[Z]=fc[Z];
                fc[Z]=fa[Z];
            }
            tol1[Z]=2.0*EPS*fabs(b[Z])+0.5*tol;
            xm[Z]=0.5*(c[Z]-b[Z]);
        }
        
        sol_check = 0;
        for(int Z =0; Z<l_Z; Z++){
            if (fabs(xm[Z]) <= tol1[Z] || fb[Z] == 0.0){sol_check += 0}else{sol_check += 1};
        }
        if(sol_check == 0){*sol1 = b1; *sol2 = b2; break;}
        
        for(int Z =0; Z<l_Z; Z++){
            if (fabs(xm[Z]) <= tol1[Z] || fb[Z] == 0.0) return b;
            if (fabs(e[Z]) >= tol1[Z] && fabs(fa[Z]) > fabs(fb[Z])) {
                s[Z]=fb[Z]/fa[Z];
                if (a[Z] == c[Z]) {
                    p[Z]=2.0*xm[Z]*s[Z];
                    q[Z]=1.0-s[Z];
                } else {
                    q[Z]=fa[Z]/fc[Z];
                    r[Z]=fb[Z]/fc[Z];
                    p[Z]=s[Z]*(2.0*xm*q*(q[Z]-r[Z])-(b[Z]-a[Z])*(r[Z]-1.0));
                    q[Z]=(q[Z]-1.0)*(r[Z]-1.0)*(s[Z]-1.0);
                }
                if (p[Z] > 0.0) q[Z] = -q[Z];
                p[Z]=fabs(p[Z]);
                min1[Z]=3.0*xm[Z]*q[Z]-fabs(tol1[Z]*q[Z]);
                min2[Z]=fabs(e[Z]*q[Z]);
                if (2.0*p[Z] < (min1[Z] < min2[Z] ? min1[Z] : min2[Z])) {
                    e[Z]=d[Z];
                    d[Z]=p[Z]/q[Z];
                } else {
                    d[Z]=xm[Z];
                    e[Z]=d[Z];
                }
            } else {
                d[Z]=xm[Z];
                e[Z]=d[Z];
            }
            a[Z]=b[Z];
            fa[Z]=fb[Z];
            if (fabs(d[Z]) > tol1[Z])
                b[Z] += d[Z];
            else
                b[Z] += SIGN(tol1[Z],xm[Z]);
            fb[Z]=func(b,params);
        }

    }
//    NR::nrerror("Maximum number of iterations exceeded in zbrent");
    return 0.0;
}



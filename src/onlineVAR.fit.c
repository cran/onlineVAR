#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <R_ext/BLAS.h>


SEXP onlineVARfit(SEXP x, SEXP y, SEXP nu0, SEXP mui, SEXP wNi, SEXP Ri, 
  SEXP ri, SEXP seq, SEXP coefi, SEXP q, SEXP abstol, SEXP trace)
{
  int N = INTEGER(GET_DIM(x))[0];
  int p = INTEGER(GET_DIM(x))[1];
  int n,i,c,d,k,Ntrace=1, ptrace = 1;
  double xn[p], yn[p], lambda, sum = 0, coefnew, diff, rpall[p*p], C[p*p];
  double *xptr = REAL(x);
  double *yptr = REAL(y);
  double *nu0ptr = REAL(nu0);
  double *muiptr = REAL(mui);
  double *wNiptr = REAL(wNi);
  double *Riptr = REAL(Ri);
  double *riptr = REAL(ri);
  double *seqptr = REAL(seq);
  double *coefiptr = REAL(coefi);
  int *qptr = INTEGER(q);
  int *traceptr = INTEGER(trace);
  double *abstolptr = REAL(abstol);
  double one = 1.0;
  double zero = 0.0;

  if(traceptr[0] == 1) { 
    Ntrace = N;
    ptrace = p;
  }
  SEXP coefall = PROTECT(allocMatrix(REALSXP, Ntrace*ptrace, ptrace));
  

  SEXP mu = PROTECT(allocVector(REALSXP, p));
  SEXP wN = PROTECT(allocVector(REALSXP, 1));
  SEXP R = PROTECT(allocMatrix(REALSXP, p, p));
  SEXP r = PROTECT(allocMatrix(REALSXP, p, p));
  SEXP coef = PROTECT(allocMatrix(REALSXP, p, p));
  SEXP pred = PROTECT(allocMatrix(REALSXP, N, qptr[0]));
  SEXP rval = PROTECT(allocVector(VECSXP, 7));
  SET_VECTOR_ELT(rval, 0, mu);
  SET_VECTOR_ELT(rval, 1, wN);
  SET_VECTOR_ELT(rval, 2, R);
  SET_VECTOR_ELT(rval, 3, r);
  SET_VECTOR_ELT(rval, 4, coef);
  SET_VECTOR_ELT(rval, 5, pred);
  SET_VECTOR_ELT(rval, 6, coefall);

  double *Rptr = REAL(R);
  double *wNptr = REAL(wN);
  double *muptr = REAL(mu);
  double *rptr = REAL(r);
  double *coefptr = REAL(coef);
  double *predptr = REAL(pred);
  double *coefallptr = REAL(coefall);

  for (i = 0; i<p; i++) {
    muptr[i] = muiptr[i];
  }
  wNptr[0] = wNiptr[0];

  for (i = 0; i<p*p; i++) {
    Rptr[i] = Riptr[i];
    rptr[i] = riptr[i];
    coefptr[i] = coefiptr[i];
  }


  for (n = 0; n < N; n++) {
/*    mu <- nu0*mu + y[n,]*/
/*    wN <- nu0*wN + 1 */
/*    yn <- if(n>2) (y[n,] - mu/wN) else y[n,]*/
/*    xn <- x[n,]- mu/wN*/
    wNptr[0] = nu0ptr[0]* wNptr[0] + 1;
    for (i = 0; i < p; i++) {
      muptr[i] = nu0ptr[0]*muptr[i] + yptr[n+N*i];
      yn[i] = yptr[n+N*i] - muptr[i]/ wNptr[0];
      xn[i] = xptr[n+N*i] - muptr[i]/ wNptr[0];

    }

    for (c = 0; c < qptr[0]; c++) {
      for (k = 0; k < p; k++) {
        sum = sum + coefptr[c + p*k]*xn[k];
      }
      predptr[n + N*c] = sum + muptr[c]/wNptr[0];
      sum = 0;
    }


/*    ## update R and r*/
/*    R <- nu0*R + xn%*%t(xn)*/
/*    r <- nu0*r + yn%*%t(xn)*/
    for (c = 0; c < p; c++) {
      for (d = 0; d < p; d++) {
        Rptr[c + p*d] = nu0ptr[0]*Rptr[c + p*d] + xn[c]*xn[d];
        rptr[c + p*d] = nu0ptr[0]*rptr[c + p*d] + yn[c]*xn[d];
      }
    }

/*    lambda <- max(abs(r))*seq*/
    lambda = 0;
    for(i = 0; i < p*p; i++) {
      lambda = fmax(fabs(rptr[i]), lambda);
    }
    lambda = lambda*seqptr[0];



/*    Rpp <- matrix(rep(diag(R), P), ncol = P, byrow = TRUE)*/
/*    rpall <- r - coef[[frac]]%*%t(R) + coef[[frac]]*Rpp*/
    F77_CALL(dgemm)("N", "T", &p,&p,&p, &one, coefptr, &p, Rptr, &p,
        &zero, C, &p);

    
    


    for (c = 0; c < p; c++) {
      for (d = 0; d < p; d++) {
        /*Only do multiplication coef not 0 speeds up significantly*/
        if(coefptr[c + p*d]!=0) {
          rpall[c + p*d] =  rptr[c + p*d] - C[c+p*d] + coefptr[c + p*d]*Rptr[d+p*d];
        } else {
          rpall[c + p*d] = rptr[c + p*d] - C[c+p*d];
        }
      }
    }

/*    for (c = 0; c < p; c++) {*/
/*      for (d = 0; d < p; d++) {*/
/*        rpall[c + p*d] = rptr[c + p*d] - C[c+p*d] + coefptr[c + p*d]*Rptr[d+p*d];*/
/*      }*/
/*    }*/


/*    for (c = 0; c < p; c++) {*/
/*      for (d = 0; d < p; d++) {*/
/*        for (k = 0; k < p; k++) {*/
/*          sum = sum + coefptr[c + p*k]*Rptr[p*k + d];*/
/*        }*/
/*        rpall[c + p*d] = rptr[c + p*d] - sum + coefptr[c + p*d]*Rptr[d+p*d];*/
/*        sum = 0;*/
/*      }*/
/*    }*/


/*     if(lambda[frac] < abs(rpall[c,k])){*/
/*       coefnew <- if(rpall[c,k]>0) (rpall[c,k]-lambda[frac])/R[k,k] else */
/*         (rpall[c,k] + lambda[frac])/R[k,k]*/
/*       rpall[c,-k] <- rpall[c,-k] - R[-k,k]*(coefnew - coef[[frac]][c,k])*/
/*       coef[[frac]][c, k] <- coefnew*/
/*     } else {*/
/*       if(coef[[frac]][c,k] != 0) {*/
/*         rpall[c,-k] <- rpall[c,-k] - R[-k,k]*(- coef[[frac]][c,k])*/
/*         coef[[frac]][c,k] <- 0*/
/*       }*/
/*     }*/

    
    diff = abstolptr[0]*p*p + 1;
    while(diff > abstolptr[0]*p*p) {

      diff = 0;
      for (c = 0; c < qptr[0]; c++) {
        for (k = 0; k < p; k++) {
          if(lambda < fabs(rpall[c + p*k])){
            if(rpall[c+p*k]>0){
              coefnew = (rpall[c+p*k]-lambda)/Rptr[k+p*k];
            } else {
              coefnew = (rpall[c+p*k]+lambda)/Rptr[k+p*k];
            }
            for(d = 0; d<p; d++) {
              if(d!=k) {
                rpall[c+p*d] = rpall[c+p*d] - Rptr[d + p*k]*(coefnew-coefptr[c+p*k]);
              }
            }
            diff = diff + fabs(coefnew-coefptr[c+p*k]);
            coefptr[c+p*k] = coefnew;
          } else {
            if(coefptr[c+p*k] != 0){
              for(d = 0; d<p; d++) {
                if(d!=k) {
                  rpall[c+p*d] = rpall[c+p*d] - Rptr[d + p*k]*(-coefptr[c+p*k]);
                }
              }
              diff = diff + fabs(coefptr[c+p*k]);
              coefptr[c+p*k] = 0;
            }
          }     
        }
      }
    }

    if(traceptr[0] == 1) {  
      for (c = 0; c < p; c++) {
        for (d = 0; d < p; d++) {
          coefallptr[c+n*p + p*N*d] = coefptr[c + p*d];
        }
      }
    }
  }




  UNPROTECT(8);
  return rval;
}







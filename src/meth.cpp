#include <R.h>
#include <Rmath.h>
#include <vector>
#include "meth.h"

void updateLambdaA ( state &Q )
{
    double newL, oldL;
    double qOldL, qNewL;
    bool jump;
    int ii;
    double t;

    for ( int i = 0; i < Q.nA; i++ )
    {
        for ( int j = 0; j < Q.m; j++ )
        {
            ii = j*Q.nA + i;
            oldL = Q.lambdaA.El[ii];
            if ( Q.step < Q.t0 )
            {
                newL = rnorm(oldL,0.1);
            }
            else
            {
                newL = rnorm(oldL,sqrt(Q.AMHLambdaA.Sigma[ii]));
            }
    
            qOldL = dnorm(oldL,Q.muA[j],sqrt(exp(Q.alphaA[j])),0)*dpois(Q.kA.El[ii],Q.sA[i]*exp(oldL),0);
            qNewL = dnorm(newL,Q.muA[j],sqrt(exp(Q.alphaA[j])),0)*dpois(Q.kA.El[ii],Q.sA[i]*exp(newL),0);

            jump = (runif(0,1) < qNewL / qOldL);
            if ( jump )
                Q.lambdaA.El[ii] = newL;
            else
                Q.lambdaA.El[ii] = oldL;
        }
    }
    for ( int i = 0; i < Q.nA; i++ )
    {
        for ( int j = 0; j < Q.m; j++ )
        {
            ii = j*Q.nA + i;
   
            // update empirical variance
            t = Q.step;
            Q.AMHLambdaA.X[ii] = Q.lambdaA.El[ii];
            Q.AMHLambdaA.Xmbar[ii] = Q.AMHLambdaA.Xbar[ii];
            Q.AMHLambdaA.Xbar[ii] = (t*Q.AMHLambdaA.Xbar[ii] + Q.lambdaA.El[ii])/(t+1);
            if ( t != 0 )
            {
                Q.AMHLambdaA.Sigma[ii] = (t-1)/t*Q.AMHLambdaA.Sigma[ii] + 2.4*2.4/t*(
                    t*Q.AMHLambdaA.Xmbar[ii]*Q.AMHLambdaA.Xmbar[ii]
                    -(t+1)*Q.AMHLambdaA.Xbar[ii]*Q.AMHLambdaA.Xbar[ii]
                    +Q.AMHLambdaA.X[ii]*Q.AMHLambdaA.X[ii]
                );
            }
        }
    }
}

void updateLambdaB ( state &Q )
{
    double newL, oldL;
    double qOldL, qNewL;
    bool jump;
    int ii;
    double t;

    for ( int i = 0; i < Q.nB; i++ )
    {
        for ( int j = 0; j < Q.m; j++ )
        {
            ii = j*Q.nB + i;
            oldL = Q.lambdaB.El[ii];
            if ( Q.step < Q.t0 )
            {
                newL = rnorm(oldL,0.1);
            }
            else
            {
                newL = rnorm(oldL,sqrt(Q.AMHLambdaB.Sigma[ii]));
            }
            qOldL = dnorm(oldL,Q.muA[j] + Q.gamma[j],sqrt(exp(Q.alphaB[j])),0)*dpois(Q.kB.El[ii],Q.sB[i]*exp(oldL),0);
            qNewL = dnorm(newL,Q.muA[j] + Q.gamma[j],sqrt(exp(Q.alphaB[j])),0)*dpois(Q.kB.El[ii],Q.sB[i]*exp(newL),0);

            jump = (runif(0,1) < qNewL / qOldL);
            if ( jump )
                Q.lambdaB.El[ii] = newL;
            else
                Q.lambdaB.El[ii] = oldL;
        }
    }
    for ( int i = 0; i < Q.nB; i++ )
    {
        for ( int j = 0; j < Q.m; j++ )
        {
            ii = j*Q.nB + i;
   
            // update empirical variance
            t = Q.step;
            Q.AMHLambdaB.X[ii] = Q.lambdaB.El[ii];
            Q.AMHLambdaB.Xmbar[ii] = Q.AMHLambdaB.Xbar[ii];
            Q.AMHLambdaB.Xbar[ii] = (t*Q.AMHLambdaB.Xbar[ii] + Q.lambdaB.El[ii])/(t+1);
            if ( t != 0 )
            {
                Q.AMHLambdaB.Sigma[ii] = (t-1)/t*Q.AMHLambdaB.Sigma[ii] + 2.4*2.4/t*(
                    t*Q.AMHLambdaB.Xmbar[ii]*Q.AMHLambdaB.Xmbar[ii]
                    -(t+1)*Q.AMHLambdaB.Xbar[ii]*Q.AMHLambdaB.Xbar[ii]
                    +Q.AMHLambdaB.X[ii]*Q.AMHLambdaB.X[ii]
                );
            }
        }
    }
}

void updateAlphaA ( state &Q )
{
    double newA, oldA;
    double qOldA, qNewA;
    bool jump;
    double t;
    int ii;

    for ( int j = 0; j < Q.m; j++ )
    {
        oldA = Q.alphaA[j];
        if ( Q.step < Q.t0 )
        {
            newA = rnorm(oldA,0.1);
        }
        else
        {
            newA = rnorm(oldA,sqrt(Q.AMHAlphaA.Sigma[j]));
        }
        qOldA = dnorm(oldA, Q.psi0, Q.tau,0);
        qNewA = dnorm(newA, Q.psi0, Q.tau,0);
        for ( int i = 0; i < Q.nA; i++ )
        {
            ii = j*Q.nA + i;
            qNewA = qNewA * dnorm(Q.lambdaA.El[ii],Q.muA[j],sqrt(exp(newA)),0);
            qOldA = qOldA * dnorm(Q.lambdaA.El[ii],Q.muA[j],sqrt(exp(oldA)),0);
        }
        jump = (runif(0,1) < qNewA / qOldA);
        if ( jump )
            Q.alphaA[j] = newA;
        else
            Q.alphaA[j] = oldA;
    }
    for ( int j = 0; j < Q.m; j++ )
    {
        // update empirical variance
        t = Q.step;
        Q.AMHAlphaA.X[j] = Q.alphaA[j];
        Q.AMHAlphaA.Xmbar[j] = Q.AMHAlphaA.Xbar[j];
        Q.AMHAlphaA.Xbar[j] = (t*Q.AMHAlphaA.Xbar[j] + Q.alphaA[j])/(t+1);
        if ( t != 0 )
        {
            Q.AMHAlphaA.Sigma[j] = (t-1)/t*Q.AMHAlphaA.Sigma[j] + 2.4*2.4/t*(
                t*Q.AMHAlphaA.Xmbar[j]*Q.AMHAlphaA.Xmbar[j]
                -(t+1)*Q.AMHAlphaA.Xbar[j]*Q.AMHAlphaA.Xbar[j]
                +Q.AMHAlphaA.X[j]*Q.AMHAlphaA.X[j]
            );
        }
    }
}

void updateAlphaB ( state &Q )
{
    double newA, oldA;
    double qOldA, qNewA;
    bool jump;
    double t;
    int ii;

    for ( int j = 0; j < Q.m; j++ )
    {
        oldA = Q.alphaB[j];
        if ( Q.step < Q.t0 )
        {
            newA = rnorm(oldA,0.1);
        }
        else
        {
            newA = rnorm(oldA,sqrt(Q.AMHAlphaB.Sigma[j]));
        }
        qOldA = dnorm(oldA, Q.psi0, Q.tau,0);
        qNewA = dnorm(newA, Q.psi0, Q.tau,0);
        for ( int i = 0; i < Q.nB; i++ )
        {
            ii = j*Q.nB + i;
            qNewA = qNewA * dnorm(Q.lambdaB.El[ii],Q.muA[j] + Q.gamma[j],sqrt(exp(newA)),0);
            qOldA = qOldA * dnorm(Q.lambdaB.El[ii],Q.muA[j] + Q.gamma[j],sqrt(exp(oldA)),0);
        }
        jump = (runif(0,1) < qNewA / qOldA);
        if ( jump )
            Q.alphaB[j] = newA;
        else
            Q.alphaB[j] = oldA;
    }
    for ( int j = 0; j < Q.m; j++ )
    {
        // update empirical variance
        t = Q.step;
        Q.AMHAlphaB.X[j] = Q.alphaB[j];
        Q.AMHAlphaB.Xmbar[j] = Q.AMHAlphaB.Xbar[j];
        Q.AMHAlphaB.Xbar[j] = (t*Q.AMHAlphaB.Xbar[j] + Q.alphaB[j])/(t+1);
        if ( t != 0 )
        {
            Q.AMHAlphaB.Sigma[j] = (t-1)/t*Q.AMHAlphaB.Sigma[j] + 2.4*2.4/t*(
                t*Q.AMHAlphaB.Xmbar[j]*Q.AMHAlphaB.Xmbar[j]
                -(t+1)*Q.AMHAlphaB.Xbar[j]*Q.AMHAlphaB.Xbar[j]
                +Q.AMHAlphaB.X[j]*Q.AMHAlphaB.X[j]
            );
        }
    }
}




void updateInd ( state &Q )
{
    double p,q;
    for ( int i = 0; i < Q.m; i++ )
    {
        p = Q.pi0 * dnorm(
                columnMean(Q.lambdaA,i),
                columnMean(Q.lambdaB,i),
               sqrt((Q.nB*exp(Q.alphaA[i]) + Q.nA*exp(Q.alphaB[i]))/(Q.nA*Q.nB) + Q.sigmaGamma*Q.sigmaGamma),0);
        q = (1-Q.pi0) * dnorm(
                columnMean(Q.lambdaA,i),
                columnMean(Q.lambdaB,i),
                sqrt((Q.nB*exp(Q.alphaA[i]) + Q.nA*exp(Q.alphaB[i]))/(Q.nA*Q.nB)),0);
        Q.ind[i] = rbinom(1,p/(p+q));
    }
}

void updateMuA ( state &Q )
{
    double mean, sd;
    for ( int i = 0; i < Q.m; i++ )
    {
        if ( Q.ind[i] == 0 )
        {
            mean = (columnMean(Q.lambdaA,i)*exp(Q.alphaB[i]) + columnMean(Q.lambdaB,i)*exp(Q.alphaA[i]))/
                    (exp(Q.alphaA[i]) + exp(Q.alphaB[i]));
            sd = sqrt(exp(Q.alphaA[i])*exp(Q.alphaB[i])/((double)Q.nB*exp(Q.alphaA[i]) + (double)Q.nA*exp(Q.alphaB[i])));
        }
        else
        {
            mean = (columnMean(Q.lambdaA,i)*(exp(Q.alphaB[i])/(double)Q.nB + Q.sigmaGamma*Q.sigmaGamma)
                         + columnMean(Q.lambdaB,i)*exp(Q.alphaA[i])/(double)Q.nA)/
                    (exp(Q.alphaA[i])/(double)Q.nA + exp(Q.alphaB[i])/(double)Q.nB + Q.sigmaGamma * Q.sigmaGamma);
            sd = sqrt(exp(Q.alphaA[i])/(double)Q.nA*(exp(Q.alphaB[i])/(double)Q.nB + Q.sigmaGamma * Q.sigmaGamma)/
                    (exp(Q.alphaA[i])/(double)Q.nA + exp(Q.alphaB[i])/(double)Q.nB + Q.sigmaGamma * Q.sigmaGamma));
        }
        Q.muA[i] = rnorm(mean,sd);
    }
}

void updateGamma ( state &Q )
{
    double mean,sd;
    for ( int i = 0; i < Q.m; i++ )
    {
        if ( Q.ind[i] == 0 )
        {
            Q.gamma[i] = 0;
        }    
        if ( Q.ind[i] != 0 )
        {
            mean = ((columnMean(Q.lambdaB,i) - Q.muA[i])*Q.sigmaGamma*Q.sigmaGamma)/
                (exp(Q.alphaB[i])/(double)Q.nB + Q.sigmaGamma*Q.sigmaGamma);
            sd = sqrt((exp(Q.alphaB[i])/(double)Q.nB * Q.sigmaGamma*Q.sigmaGamma)/
                    (exp(Q.alphaB[i])/(double)Q.nB + Q.sigmaGamma*Q.sigmaGamma));
            Q.gamma[i] = rnorm(mean,sd);
        }
    }
}

void updatePi0 ( state &Q )
{
    Q.pi0 = rbeta(sum(Q.ind) + 1, (double)Q.m - sum(Q.ind) + 1);
}

void updateSigmaGamma ( state &Q )
{
    double shape = 0.5*sum(Q.ind);
    double rate = 0.5*sum(Q.gamma*Q.gamma);
    Q.sigmaGamma = sqrt(rinvgamma(shape,rate));
}

void updatePsi0 ( state &Q )
{
    double mean = 1/(2*(double)Q.m/(Q.tau*Q.tau))*(sum(Q.alphaA)/(Q.tau*Q.tau) + sum(Q.alphaB)/(Q.tau*Q.tau));
    double sd = sqrt(1/(2*(double)Q.m/(Q.tau*Q.tau)));
    Q.psi0 = rnorm(mean,sd);
}

void updateTau ( state &Q )
{
    double alpha, beta;
    alpha = Q.m;
    beta = 0;
    for ( int i = 0; i < Q.m; i ++ )
    {
        beta += (Q.alphaA[i] - Q.psi0)*(Q.alphaA[i] - Q.psi0);
        beta += (Q.alphaB[i] - Q.psi0)*(Q.alphaB[i] - Q.psi0);
    }
    beta = 0.5 * beta;
    Q.tau = sqrt(rinvgamma(alpha,beta));
}

void updateState ( state &Q )
{
    GetRNGstate();
    updateLambdaA ( Q );
    updateLambdaB ( Q );
    updateAlphaA ( Q );
    updateAlphaB ( Q );
    updateInd ( Q );
    updateMuA ( Q );
    updateGamma ( Q );
    updatePi0 ( Q );
    updateSigmaGamma ( Q );
    updatePsi0 ( Q );
    updateTau ( Q );

    Q.step++;
    PutRNGstate();
}

extern "C" {

    void rnaseq ( double *kA, double *kB, double *sA, double *sB, int *nA, int *nB, int *m,
        int *burn, int *reps, int *saveEvery, int *printEvery, int *t0, 
        double *lambdaA, double *lambdaB, double *ind, double *muA, double *gamma, //
        double *alphaA, double *alphaB, double *pi0, double *sigmaGamma,           // start parameters
        double *psi0, double *tau, int *success )                                                //
    {

        state Q;

        // initialize state
        Q.step = 0;

        Q.nA = *nA;
        Q.nB = *nB;
        Q.m = *m;

        Q.kA.El.assign(kA,kA+*nA* *m);
        Q.kA.rows = *nA;
        Q.kA.columns = *m;

        Q.kB.El.assign(kB,kB+*nB* *m);
        Q.kB.rows = *nB;
        Q.kB.columns = *m;

        Q.sA.assign(sA, sA + *nA);
        Q.sB.assign(sB, sB + *nB);
        
        Q.lambdaA.El.assign(lambdaA,lambdaA + *nA * *m);
        Q.lambdaA.rows = *nA;
        Q.lambdaA.columns = *m;

        Q.AMHLambdaA.X.assign(lambdaA, lambdaA + *nA * *m); 
        Q.AMHLambdaA.Xbar.assign(lambdaA, lambdaA + *nA * *m);
        Q.AMHLambdaA.Xmbar.assign(*nA * *m,0);
        Q.AMHLambdaA.Sigma.assign(*nA * *m,0);

        Q.lambdaB.El.assign(lambdaB,lambdaB + *nB * *m);
        Q.lambdaB.rows = *nB;
        Q.lambdaB.columns = *m;

        Q.AMHLambdaB.X.assign(lambdaB, lambdaB + *nB * *m); 
        Q.AMHLambdaB.Xbar.assign(lambdaB, lambdaB + *nB * *m);
        Q.AMHLambdaB.Xmbar.assign(*nB * *m,0);
        Q.AMHLambdaB.Sigma.assign(*nB * *m,0);

        Q.alphaA.assign(alphaA, alphaA + *m);

        Q.AMHAlphaA.X.assign(alphaA, alphaA + *m); 
        Q.AMHAlphaA.Xbar.assign(alphaA, alphaA + *m);
        Q.AMHAlphaA.Xmbar.assign(*m,0);
        Q.AMHAlphaA.Sigma.assign(*m,0);

        Q.alphaB.assign(alphaB, alphaB + *m);
          
        Q.AMHAlphaB.X.assign(alphaB, alphaB + *m); 
        Q.AMHAlphaB.Xbar.assign(alphaB, alphaB + *m);
        Q.AMHAlphaB.Xmbar.assign(*m,0);
        Q.AMHAlphaB.Sigma.assign(*m,0);

        Q.ind.assign(ind, ind + *m);
        Q.muA.assign(muA, muA + *m);
        Q.gamma.assign(gamma, gamma + *m);
        
        Q.pi0 = *pi0;
        Q.psi0 = *psi0;
        
        Q.sigmaGamma = *sigmaGamma;
        Q.tau = *tau;

        Rprintf("++++++++++++++++++++++++++++++++++++++\n");        

        for ( int i = 0; i < *burn; i++ )
        {
            if ( (i + 1) % *printEvery == 0 )
                Rprintf("++++++++++ Burn Step %6d ++++++++++\n", i + 1);
           
            updateState(Q);

            if ( Q.sigmaGamma != Q.sigmaGamma )
                return;

            R_CheckUserInterrupt();
        }

        Rprintf("++++++++++++++++++++++++++++++++++++++\n");        

        for ( int i = 0; i < *reps; i++ )
        {
            if ( (i + 1) % *printEvery == 0 )
                Rprintf("++++++++++ Repetition %5d ++++++++++\n", i + 1 );

            updateState(Q);

            if ( Q.sigmaGamma != Q.sigmaGamma )
                return;

            R_CheckUserInterrupt();
            // save results //

            if ( i == 0 )
            {
                for ( int j = 0; j < Q.m; j++ )
                {
/*
                    for ( int k = 0; k < Q.nA; k++ )
                    {
                        ii = j*Q.nA + k;
                        *(lambdaA + ii) = Q.lambdaA.El[ii];
                    }
                    for ( int k = 0; k < Q.nB; k++ )
                    {
                        ii = j*Q.nB + k;
                        *(lambdaB + ii) = Q.lambdaB.El[ii];
                    }
*/
                    *(alphaA + j) = Q.alphaA[j];
                    *(alphaB + j) = Q.alphaB[j];
                    *(ind + j) = Q.ind[j];
                    *(gamma + j) = Q.gamma[j];
                    *(muA +j) = Q.muA[j];
                }
                *pi0 = Q.pi0;
                *sigmaGamma = Q.sigmaGamma;
                *psi0 = Q.psi0;
                *tau = Q.tau;
            }
            else
            {
                for ( int j = 0; j < Q.m; j++ )
                {
/*
                    for ( int k = 0; k < Q.nA; k++ )
                    {
                        ii = j*Q.nA + k;
                        *(lambdaA + ii) = (*(lambdaA+ii) * i + Q.lambdaA.El[ii])/(i+1);
                    }
                    for ( int k = 0; k < Q.nB; k++ )
                    {
                        ii = j*Q.nB + k;
                        *(lambdaB + ii) = (*(lambdaB+ii) * i + Q.lambdaB.El[ii])/(i+1);
                    }
*/

                    *(alphaA + j) = (*(alphaA+j) * i + Q.alphaA[j])/(i+1);
                    *(alphaB + j) = (*(alphaB+j) * i + Q.alphaB[j])/(i+1);
                    *(ind + j) = (*(ind+j) * i + Q.ind[j])/(i+1);
                    *(gamma + j) = (*(gamma+j) * i + Q.gamma[j])/(i+1);
                    *(muA + j) = (*(muA+j) * i + Q.muA[j])/(i+1);

                }
                *pi0 = (*pi0 * i + Q.pi0)/(i+1);

                *sigmaGamma = (*sigmaGamma * i + Q.sigmaGamma)/(i+1);
                *psi0 = (*psi0 * i + Q.psi0)/(i+1);
                *tau = (*tau * i + Q.tau)/(i+1);
            }

        }
        Rprintf("++++++++++++++++++++++++++++++++++++++\n");        
        *success = 1;
    }


    void rnaseq_post_dist ( double *kA, double *kB, double *sA, double *sB, int *nA, int *nB, int *m,
        int *burn, int *reps, int *saveEvery, int *printEvery, int *t0, 
        double *lambdaA, double *lambdaB, double *ind, double *muA, double *gamma, //
        double *alphaA, double *alphaB, double *pi0, double *sigmaGamma,           // start parameters
        double *psi0, double *tau,                                                 //
        double *distGamma, int *success )                                                
    {
        int ii, kk;

        state Q;

        // initialize state
        Q.step = 0;

        Q.nA = *nA;
        Q.nB = *nB;
        Q.m = *m;

        Q.kA.El.assign(kA,kA+*nA* *m);
        Q.kA.rows = *nA;
        Q.kA.columns = *m;

        Q.kB.El.assign(kB,kB+*nB* *m);
        Q.kB.rows = *nB;
        Q.kB.columns = *m;

        Q.sA.assign(sA, sA + *nA);
        Q.sB.assign(sB, sB + *nB);
        
        Q.lambdaA.El.assign(lambdaA,lambdaA + *nA * *m);
        Q.lambdaA.rows = *nA;
        Q.lambdaA.columns = *m;

        Q.AMHLambdaA.X.assign(lambdaA, lambdaA + *nA * *m); 
        Q.AMHLambdaA.Xbar.assign(lambdaA, lambdaA + *nA * *m);
        Q.AMHLambdaA.Xmbar.assign(*nA * *m,0);
        Q.AMHLambdaA.Sigma.assign(*nA * *m,0);

        Q.lambdaB.El.assign(lambdaB,lambdaB + *nB * *m);
        Q.lambdaB.rows = *nB;
        Q.lambdaB.columns = *m;

        Q.AMHLambdaB.X.assign(lambdaB, lambdaB + *nB * *m); 
        Q.AMHLambdaB.Xbar.assign(lambdaB, lambdaB + *nB * *m);
        Q.AMHLambdaB.Xmbar.assign(*nB * *m,0);
        Q.AMHLambdaB.Sigma.assign(*nB * *m,0);

        Q.alphaA.assign(alphaA, alphaA + *m);

        Q.AMHAlphaA.X.assign(alphaA, alphaA + *m); 
        Q.AMHAlphaA.Xbar.assign(alphaA, alphaA + *m);
        Q.AMHAlphaA.Xmbar.assign(*m,0);
        Q.AMHAlphaA.Sigma.assign(*m,0);

        Q.alphaB.assign(alphaB, alphaB + *m);
          
        Q.AMHAlphaB.X.assign(alphaB, alphaB + *m); 
        Q.AMHAlphaB.Xbar.assign(alphaB, alphaB + *m);
        Q.AMHAlphaB.Xmbar.assign(*m,0);
        Q.AMHAlphaB.Sigma.assign(*m,0);

        Q.ind.assign(ind, ind + *m);
        Q.muA.assign(muA, muA + *m);
        Q.gamma.assign(gamma, gamma + *m);
        
        Q.pi0 = *pi0;
        Q.psi0 = *psi0;
        
        Q.sigmaGamma = *sigmaGamma;
        Q.tau = *tau;

        Rprintf("++++++++++++++++++++++++++++++++++++++\n");        

        for ( int i = 0; i < *burn; i++ )
        {
            if ( (i + 1) % *printEvery == 0 )
                Rprintf("++++++++++ Burn Step %6d ++++++++++\n", i + 1);
           
            updateState(Q);

            if ( Q.sigmaGamma != Q.sigmaGamma )
                return;

            R_CheckUserInterrupt();
        }

        Rprintf("++++++++++++++++++++++++++++++++++++++\n");        

        kk = 0;
        for ( int i = 0; i < *reps; i++ )
        {
            if ( (i + 1) % *printEvery == 0 )
                Rprintf("++++++++++ Repetition %5d ++++++++++\n", i + 1 );

            updateState(Q);

            if ( Q.sigmaGamma != Q.sigmaGamma )
                return;

            R_CheckUserInterrupt();
            // save results //

            if ( i % *saveEvery == 0)
            {
                if ( kk == 0 )
                {
                    for ( int j = 0; j < Q.m; j++ )
                    {
    /*
                        for ( int k = 0; k < Q.nA; k++ )
                        {
                            ii = j*Q.nA + k;
                            *(lambdaA + ii) = Q.lambdaA.El[ii];
                        }
                        for ( int k = 0; k < Q.nB; k++ )
                        {
                            ii = j*Q.nB + k;
                            *(lambdaB + ii) = Q.lambdaB.El[ii];
                        }
    */
                        *(alphaA + j) = Q.alphaA[j];
                        *(alphaB + j) = Q.alphaB[j];
                        *(ind + j) = Q.ind[j];
                        *(gamma + j) = Q.gamma[j];
                        *(muA +j) = Q.muA[j];
                    }
                    *pi0 = Q.pi0;
                    *sigmaGamma = Q.sigmaGamma;
                    *psi0 = Q.psi0;
                    *tau = Q.tau;
                }
                else {
                    for ( int j = 0; j < Q.m; j++ )
                    {
  /*
                        for ( int k = 0; k < Q.nA; k++ )
                        {
                            ii = j*Q.nA + k;
                            *(lambdaA + ii) = (*(lambdaA+ii) * i + Q.lambdaA.El[ii])/(i+1);
                        }
                        for ( int k = 0; k < Q.nB; k++ )
                        {
                            ii = j*Q.nB + k;
                            *(lambdaB + ii) = (*(lambdaB+ii) * i + Q.lambdaB.El[ii])/(i+1);
                        }
    */

                        *(alphaA + j) = (*(alphaA+j) * kk + Q.alphaA[j])/(kk+1);
                        *(alphaB + j) = (*(alphaB+j) * kk + Q.alphaB[j])/(kk+1);
                        *(ind + j) = (*(ind+j) * kk + Q.ind[j])/(kk+1);
                        *(gamma + j) = (*(gamma+j) * kk + Q.gamma[j])/(kk+1);
                        *(muA + j) = (*(muA+j) * kk + Q.muA[j])/(kk+1);

                    }
                    *pi0 = (*pi0 * kk + Q.pi0)/(kk+1);

                    *sigmaGamma = (*sigmaGamma * kk + Q.sigmaGamma)/(kk+1);
                    *psi0 = (*psi0 * kk + Q.psi0)/(kk+1);
                    *tau = (*tau * kk + Q.tau)/(kk+1);
                }
                // save posterior distribution of gamma
                for ( int j = 0; j < Q.m; j++ )
                {
                    ii = j*(*reps / *saveEvery) + kk;
                    *(distGamma + ii) = Q.gamma[j];
                }
                kk = kk + 1;

            }
        }
        Rprintf("++++++++++++++++++++++++++++++++++++++\n");        

        *success = 1;
    }

    void rnaseq_verbose ( double *kA, double *kB, double *sA, double *sB, int *nA, int *nB, int *m,
        int *burn, int *reps, int *saveEvery, int *printEvery, int *t0, 
        double *lambdaA, double *lambdaB, double *ind, double *muA, double *gamma, //
        double *alphaA, double *alphaB, double *pi0, double *sigmaGamma,           // start parameters
        double *psi0, double *tau,                                                 //
        double *distGamma, double *distMuA, double *distAlphaA, double *distAlphaB,
        int *success
    )                  
    {
        int ii, kk;

        state Q;

        // initialize state
        Q.step = 0;

        Q.nA = *nA;
        Q.nB = *nB;
        Q.m = *m;

        Q.kA.El.assign(kA,kA+*nA* *m);
        Q.kA.rows = *nA;
        Q.kA.columns = *m;

        Q.kB.El.assign(kB,kB+*nB* *m);
        Q.kB.rows = *nB;
        Q.kB.columns = *m;

        Q.sA.assign(sA, sA + *nA);
        Q.sB.assign(sB, sB + *nB);
        
        Q.lambdaA.El.assign(lambdaA,lambdaA + *nA * *m);
        Q.lambdaA.rows = *nA;
        Q.lambdaA.columns = *m;

        Q.AMHLambdaA.X.assign(lambdaA, lambdaA + *nA * *m); 
        Q.AMHLambdaA.Xbar.assign(lambdaA, lambdaA + *nA * *m);
        Q.AMHLambdaA.Xmbar.assign(*nA * *m,0);
        Q.AMHLambdaA.Sigma.assign(*nA * *m,0);

        Q.lambdaB.El.assign(lambdaB,lambdaB + *nB * *m);
        Q.lambdaB.rows = *nB;
        Q.lambdaB.columns = *m;

        Q.AMHLambdaB.X.assign(lambdaB, lambdaB + *nB * *m); 
        Q.AMHLambdaB.Xbar.assign(lambdaB, lambdaB + *nB * *m);
        Q.AMHLambdaB.Xmbar.assign(*nB * *m,0);
        Q.AMHLambdaB.Sigma.assign(*nB * *m,0);

        Q.alphaA.assign(alphaA, alphaA + *m);

        Q.AMHAlphaA.X.assign(alphaA, alphaA + *m); 
        Q.AMHAlphaA.Xbar.assign(alphaA, alphaA + *m);
        Q.AMHAlphaA.Xmbar.assign(*m,0);
        Q.AMHAlphaA.Sigma.assign(*m,0);

        Q.alphaB.assign(alphaB, alphaB + *m);
          
        Q.AMHAlphaB.X.assign(alphaB, alphaB + *m); 
        Q.AMHAlphaB.Xbar.assign(alphaB, alphaB + *m);
        Q.AMHAlphaB.Xmbar.assign(*m,0);
        Q.AMHAlphaB.Sigma.assign(*m,0);

        Q.ind.assign(ind, ind + *m);
        Q.muA.assign(muA, muA + *m);
        Q.gamma.assign(gamma, gamma + *m);
        
        Q.pi0 = *pi0;
        Q.psi0 = *psi0;
        
        Q.sigmaGamma = *sigmaGamma;
        Q.tau = *tau;

        Rprintf("++++++++++++++++++++++++++++++++++++++\n");        

        for ( int i = 0; i < *burn; i++ )
        {
            if ( (i + 1) % *printEvery == 0 )
                Rprintf("++++++++++ Burn Step %6d ++++++++++\n", i + 1);
           
            updateState(Q);

            if ( Q.sigmaGamma != Q.sigmaGamma )
                return;

            R_CheckUserInterrupt();
        }

        Rprintf("++++++++++++++++++++++++++++++++++++++\n");        

        kk = 0;
        for ( int i = 0; i < *reps; i++ )
        {
            if ( (i + 1) % *printEvery == 0 )
                Rprintf("++++++++++ Repetition %5d ++++++++++\n", i + 1 );

            updateState(Q);

            if ( Q.sigmaGamma != Q.sigmaGamma )
                return;

            R_CheckUserInterrupt();
            // save results //

            if ( i % *saveEvery == 0)
            {
                if ( kk == 0 )
                {
                    for ( int j = 0; j < Q.m; j++ )
                    {
    /*
                        for ( int k = 0; k < Q.nA; k++ )
                        {
                            ii = j*Q.nA + k;
                            *(lambdaA + ii) = Q.lambdaA.El[ii];
                        }
                        for ( int k = 0; k < Q.nB; k++ )
                        {
                            ii = j*Q.nB + k;
                            *(lambdaB + ii) = Q.lambdaB.El[ii];
                        }
    */
                        *(alphaA + j) = Q.alphaA[j];
                        *(alphaB + j) = Q.alphaB[j];
                        *(ind + j) = Q.ind[j];
                        *(gamma + j) = Q.gamma[j];
                        *(muA +j) = Q.muA[j];
                    }
                    *pi0 = Q.pi0;
                    *sigmaGamma = Q.sigmaGamma;
                    *psi0 = Q.psi0;
                    *tau = Q.tau;
                }
                else {
                    for ( int j = 0; j < Q.m; j++ )
                    {
  /*
                        for ( int k = 0; k < Q.nA; k++ )
                        {
                            ii = j*Q.nA + k;
                            *(lambdaA + ii) = (*(lambdaA+ii) * i + Q.lambdaA.El[ii])/(i+1);
                        }
                        for ( int k = 0; k < Q.nB; k++ )
                        {
                            ii = j*Q.nB + k;
                            *(lambdaB + ii) = (*(lambdaB+ii) * i + Q.lambdaB.El[ii])/(i+1);
                        }
    */

                        *(alphaA + j) = (*(alphaA+j) * kk + Q.alphaA[j])/(kk+1);
                        *(alphaB + j) = (*(alphaB+j) * kk + Q.alphaB[j])/(kk+1);
                        *(ind + j) = (*(ind+j) * kk + Q.ind[j])/(kk+1);
                        *(gamma + j) = (*(gamma+j) * kk + Q.gamma[j])/(kk+1);
                        *(muA + j) = (*(muA+j) * kk + Q.muA[j])/(kk+1);

                    }
                    *pi0 = (*pi0 * kk + Q.pi0)/(kk+1);

                    *sigmaGamma = (*sigmaGamma * kk + Q.sigmaGamma)/(kk+1);
                    *psi0 = (*psi0 * kk + Q.psi0)/(kk+1);
                    *tau = (*tau * kk + Q.tau)/(kk+1);
                }
                // save posterior distributions
                for ( int j = 0; j < Q.m; j++ )
                {
                    ii = j*(*reps / *saveEvery) + kk;
                    *(distGamma + ii) = Q.gamma[j];
                    *(distMuA + ii) = Q.muA[j];
                    *(distAlphaA + ii ) = Q.alphaA[j];
                    *(distAlphaB + ii ) = Q.alphaB[j];
                }

                kk = kk + 1;

            }
        }
        Rprintf("++++++++++++++++++++++++++++++++++++++\n");        
        *success = 1;

    }


}

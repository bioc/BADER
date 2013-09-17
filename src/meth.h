typedef std::vector<double> Vector;
typedef struct MatrixStructure {
    std::vector<double> El;
    int rows;
    int columns;
} Matrix;

typedef struct AMHStructure {
    Vector X;
    Vector Xbar;
    Vector Xmbar;
    Vector Sigma;
} AMH;

typedef struct StateStructure {
    int step;
    int nA;
    int nB;
    int m;
    Matrix kA;
    Matrix kB;
    Vector sA;
    Vector sB;
    Matrix lambdaA;
    AMH AMHLambdaA;
    Matrix lambdaB;
    AMH AMHLambdaB;
    Vector ind;
    Vector muA;
    Vector gamma;
    Vector alphaA;
    AMH AMHAlphaA;
    Vector alphaB;
    AMH AMHAlphaB;
    double pi0;
    double sigmaGamma;
    double psi0;
    double tau;
    int t0;
} state;


// column is assumed to be from 0 to n - 1 
double columnMean ( Matrix &A, int column )
{
    double s = 0;
    for ( int i = 0; i < A.rows; i ++ )
    {
        s += A.El[column*A.rows + i];
    }
    return s / A.rows;
}

double sum ( Vector X )
{
    double s=0;
    for ( int i = 0; i < (int)X.size(); i ++ )
    {
        s += X[i];
    }
    return s;
}

double rinvgamma ( double shape, double rate )
{
    return 1/rgamma(shape,1/rate);
}

void rnorm ( Vector &mu, Vector &sd, Vector &X )
{
    for ( int i = 0; i < (int)X.size(); i ++ )
    {
        X[i] = rnorm(mu[i],sd[i]);
    }
}

Vector rnorm ( Vector &mu, Vector &sd)
{
    Vector X = mu;
    for ( int i = 0; i < (int)X.size(); i ++ )
    {
        X[i] += rnorm(0,sd[i]);
    }
    return X;
}

void rgamma ( Vector &alpha, Vector &beta, Vector &X )
{
    for ( int i = 0; i < (int)X.size(); i ++ )
    {
        X[i] = rgamma(alpha[i],beta[i]);
    }
}

void mult ( Vector &a, Matrix &X, Vector &b )
{
    int n = X.rows;
    int m = X.columns;

    for ( int j = 0; j < m; j ++ )
    {
        b[j] = 0;
        for ( int i = 0; i < n; i ++ )
        {
            b[j] += a[i] * X.El[n*j+i];
        }
    }
}

Vector operator* ( Vector &a, Vector &b )
{
    Vector c = a;
    for ( int i = 0; i < (int)a.size(); i ++ )
    {
        c[i] = a[i] * b[i];
    }
    return c;
}

Vector operator+ ( Vector &a, Vector &b )
{
    Vector c = a;
    for ( int i = 0; i < (int)a.size(); i ++ )
    {
        c[i] = a[i] + b[i];
    }
    return c;
}



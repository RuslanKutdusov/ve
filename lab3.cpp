#include <iostream>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <functional>

typedef double FP;
const FP e = 2.718281828;

void sweep( FP** a, FP* b, int N )
{
	FP znam; 

	b[ 0 ] /= a[ 0 ][ 0 ];//Q1
	a[ 0 ][ 1 ] /= -a[ 0 ][ 0 ];//P1

	for( int i = 1; i < N-1; i++ )
	{
		znam = -a[ i ][ i ] - a[ i ][ i - 1 ] * a[ i - 1 ][ i ]; //общий знаменатель для формул нахождения Pi, Qi 
		a[ i ][ i + 1 ] /= znam; //Pi
		b[ i ] = ( a[ i ][ i - 1 ] * b[ i - 1 ] - b[ i ] ) / znam; //Qi
	}
	//строка ниже для вычисления QN
	b[ N - 1 ] = ( a [ N - 1 ][ N - 2 ] * b[ N - 2 ] - b[ N - 1 ]) / ( -a[ N - 1][ N - 1 ] - a[ N - 1 ][ N - 2 ] * a[ N - 2 ][ N - 1 ]);

	//обратный ход
	for( int i = N - 2; i > -1; i-- )
	{
		b[ i ] += b[ i + 1 ] * a[ i ][ i + 1 ];
	}

	return;
}


//
void TridiagonalMatrixAlgorithm( int N, FP* solution )
{
	FP h = 1.0 / ( FP )N;

	auto di = [h]( FP xi )
	{
		FP h2 = h * h;
		return -2.9 * xi * xi * h2 + 2.9 * xi * h2 + 7.8 * h2;
	};

	auto dn = [h]()
	{
		return -e * h + h / e - 2.9 * h;
	};

	FP** A = new FP*[ N ];

	FP* B = solution;

	FP xi = 0.0;

	for( int i = 0; i < N - 1; i++ )
	{
		A[ i ] = new FP[ N ];

		FP b = 2.0 + h * h;

		if( i != 0 )
			A[ i ][ i - 1 ] = 1.0;

		A[ i ][ i ] 	= -b;
		A[ i ][ i + 1 ] = 1.0;

		B[ i ] = di( xi );
		xi += h;
	}

	A[ N - 1 ] = new FP[ N ];
	
	A[ N - 1 ][ N - 2 ] = 1.0;
	A[ N - 1 ][ N - 1 ] = -1.0;
	B[ N - 1 ] = dn();

	sweep( A, B, N );
}


//
typedef std::function< FP( FP, FP, FP ) > Func;


//
double RungeKutta( Func& f, FP h, FP xi, FP yi, FP zi )
{
	FP k1 = f( xi, 				yi, 			zi );
	FP k2 = f( xi + h / 2.0, 	yi + k1 / 2.0, 	zi + k1 / 2.0 );
	FP k3 = f( xi + h / 2.0,	yi + k2 / 2.0, 	zi + k2 / 2.0 );
	FP k4 = f( xi + h, 			yi + k3, 		zi + k3 );

	return h * ( k1 + 2.0 * k2 + 2.0 * k3 + k4 ) / 6.0;	
}


//
void RungeKutta( Func& V, Func& U, FP h, int N, FP x0, FP* yi, FP* zi )
{
	FP xi = x0;

	for( int i = 1; i < N; i++ )
	{
		yi[ i ] = yi[ i - 1 ] + RungeKutta( V, h, xi, yi[ i - 1 ], zi[ i - 1 ] );
		zi[ i ] = zi[ i - 1 ] + RungeKutta( U, h, xi, yi[ i - 1 ], zi[ i - 1 ] );

		xi += h;
	}	
}


//
void Shoot( const int N, FP* yiRunge )
{
	FP x0 = 0.0;
	FP h = 1.0 / ( FP )N;
	FP y0 = 0.0;
	FP z1 = e - 1.0 / e + 2.9;

	FP ziRunge[ N ];
	FP A[ 2 ] = { -10.0, 10.0 };
	FP B[ 2 ];

	Func V = []( FP x, FP y, FP z ) -> FP
	{
		return z;
	};

	Func U = []( FP x, FP y, FP z ) -> FP
	{
		return y + 7.8 - 2.9 * x * x + 2.9 * x;
	};

	yiRunge[ 0 ] = y0;

	ziRunge[ 0 ] = A[ 0 ];
	RungeKutta( V, U, h, N, x0, yiRunge, ziRunge );

	B[ 0 ] = ziRunge[ N - 1 ];

	while( 1 )
	{
		ziRunge[ 0 ] = A[ 1 ];
		RungeKutta( V, U, h, N, x0, yiRunge, ziRunge );

		B[ 1 ] = ziRunge[ N - 1 ];
		printf( "%f %f %f %f\n", A[ 0 ], A[ 1 ], B[ 0 ], B[ 1 ] );
		if( fabs( B[ 1 ] - z1 ) < 0.0001 || isnan( B[ 1 ] ) )
		{
			printf("%f\n", B[ 1 ] );
			break;
		}

		FP A2 = A[ 1 ] - ( B[ 1 ] - z1 ) * ( A[ 1 ] - A[ 0 ] ) / ( ( B[ 1 ] - z1 ) - ( B[ 0 ] - z1 ) );

		A[ 0 ] = A[ 1 ];
		A[ 1 ] = A2;
		B[ 0 ] = B[ 1 ];
	}
}


//
int main( int argc, char* argv[] )
{
	if( argc != 2 )
		return -1;

	int N = atoi( argv[ 1 ] );
	FP h = 1.0 / ( FP )N;

	FP* solution = new FP[ N ];

	TridiagonalMatrixAlgorithm( N, solution );

	FILE* f = fopen( "data", "w" );

	FP xi = 0.0f;

	for( int i = 0; i < N; i++ )
	{
		FP y = -2.0 + exp( -xi ) + exp( xi ) - 2.9 * xi + 2.9 * xi * xi;
		fprintf( f, "%.06f %.06f %.06f\n", xi, y, solution[ i ] );

		xi += h;
	}

	return 0;
}
#include <iostream>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <functional>

#include <vector>
#include <memory>


typedef double FP;
const FP e = 2.718281828;

void sweep( std::vector< std::vector< FP > >& a, std::vector< FP >& b, int N )
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
void TridiagonalMatrixAlgorithm( int N, std::vector< FP >& solution )
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

	std::vector< std::vector< FP > > A( N + 1, std::vector< FP >( N + 1 ) );
	auto& B = solution;

	FP xi = 0.0;

	for( int i = 0; i < N; i++ )
	{
		FP b = 2.0 + h * h;

		if( i != 0 )
			A[ i ][ i - 1 ] = 1.0;

		A[ i ][ i ] 	= -b;
		A[ i ][ i + 1 ] = 1.0;

		B[ i ] = di( xi );
		xi += h;
	}
	
	A[ N ][ N - 1 ] = 1.0;
	A[ N ][ N ] = -1.0;
	B[ N ] = dn();

	sweep( A, B, N + 1 );
}


//
void TridiagonalMatrixAlgorithmLagrange( int N, std::vector< FP >& solution )
{
	FP h = 1.0 / ( FP )N;

	auto di = [h]( FP xi )
	{
		FP h2 = h * h;
		return -2.9 * xi * xi * h2 + 2.9 * xi * h2 + 7.8 * h2;
	};

	auto dn = [h]()
	{
		return 2.0 * h * ( e - 1.0 / e + 2.9 );
	};

	std::vector< std::vector< FP > > A( N + 1, std::vector< FP >( N + 1 ) );
	auto& B = solution;

	FP xi = 0.0;

	FP b = 2.0 + h * h;

	for( int i = 0; i < N; i++ )
	{
		if( i != 0 )
			A[ i ][ i - 1 ] = 1.0;

		A[ i ][ i ] 	= -b;
		A[ i ][ i + 1 ] = 1.0;

		B[ i ] = di( xi );
		xi += h;
	}
	
	A[ N ][ N - 1 ] = b - 4.0;
	A[ N ][ N ] = 2.0;
	B[ N ] = dn() - di( 1.0 - h );

	sweep( A, B, N + 1 );
}


//
typedef std::function< FP( FP, FP, FP ) > Func;


//
void RungeKutta( Func& V, Func& U, FP h, int N, FP x0, FP* y, FP* z )
{
	FP xi = x0;

	for( int i = 0; i < N; i++ )
	{
		const FP& zi = z[ i ];
		const FP& yi = y[ i ];

		FP q0 = U( xi, 				yi,		 			zi );
		FP k0 = V( 0.0, 			0.0, 				zi );

		FP q1 = U( xi + h / 2.0, 	yi + k0 * h / 2.0, 	zi + q0 * h / 2.0 );
		FP k1 = V( 0.0, 			0.0, 				zi + q0 * h / 2.0 );

		FP q2 = U( xi + h / 2.0,	yi + k1 * h / 2.0, 	zi + q1 * h / 2.0 );
		FP k2 = V( 0.0, 			0.0, 				zi + q1 * h / 2.0 );

		FP q3 = U( xi + h,			yi + k2 * h,	 	zi + q2 * h );
		FP k3 = V( 0.0, 			0.0, 				zi + q2 * h );

		y[ i + 1 ] =  yi + h * ( k0 + 2.0 * k1 + 2.0 * k2 + k3 ) / 6.0;	
		z[ i + 1 ] =  zi + h * ( q0 + 2.0 * q1 + 2.0 * q2 + q3 ) / 6.0;	

		xi += h;
	}	
}


//
void Shoot( const int N, std::vector< FP >& yiRunge )
{
	FP x0 = 0.0;
	FP h = 1.0 / ( FP )N;
	FP y0 = 0.0;
	FP z1 = e - 1.0 / e + 2.9;

	Func V = []( FP x, FP y, FP z ) -> FP
	{
		return z;
	};

	Func U = []( FP x, FP y, FP z ) -> FP
	{
		return y + 7.8 - 2.9 * x * x + 2.9 * x;
	};

	// Рунге-Кутта
	FP eta = 0.0;

	yiRunge[ 0 ] = y0;
	FP ziRunge[ N + 1 ];

	while( 1 )
	{
		ziRunge[ 0 ] = eta;

		RungeKutta( V, U, h, N, x0, yiRunge.data(), ziRunge );

		if( fabs( ziRunge[ N ] - z1 ) < 0.1 ) 
				break;

		eta -= ( ziRunge[ N ] - z1 ) / ziRunge[ N ];
	}

	// Эйлер
}


//
int main( int argc, char* argv[] )
{
	if( argc != 2 )
		return -1;

	int N = atoi( argv[ 1 ] );
	FP h = 1.0 / ( FP )N;

	std::vector< FP > solution1( N + 1 );
	std::vector< FP > solution2( N + 1 );
	std::vector< FP > solution3( N + 1 );

	TridiagonalMatrixAlgorithm( N, solution1 );
	TridiagonalMatrixAlgorithmLagrange( N, solution2 );
	Shoot( N, solution3 );

	FILE* f = fopen( "data", "w" );

	FP xi = 0.0f;

	for( int i = 0; i <= N; i++ )
	{
		FP y = -2.0 + exp( -xi ) + exp( xi ) - 2.9 * xi + 2.9 * xi * xi;
		fprintf( f, "%.06f %.06f %.06f %.06f %.06f\n", xi, y, solution1[ i ], solution2[ i ], solution3[ i ] );

		xi += h;
	}

	return 0;
}
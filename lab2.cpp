#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <functional>

// 0 Эйлер
// 1 Эйлер с пересчетом
// 2 Эйлер неявный
// 3 Рунге-Кутта
// 4 Коши
// 5 Адамс 2 шаговый
// 6 Адамс 3 шаговый
// 7 Симпсона

std::vector< double > methods[ 8 ];
int n;
double g_y0 = 0.1;
double h;

double f( double x, double y )
{
	return 50.0 * y * ( x - 0.6 ) * ( x - 0.85 );
}


// Эйлер
void Euler()
{
	double yi = g_y0;
	double xi = 0.0f;

	methods[ 0 ].push_back( g_y0 );

	for( int i = 1; i < n; i++ )
	{
		yi = yi + h * f( xi, yi );
		methods[ 0 ].push_back( yi );

		xi += h;
	}
}


// Эйлер с пересчетом
void Euler2W()
{
	double yi = g_y0;
	double xi = 0.0f;

	methods[ 1 ].push_back( g_y0 );

	for( int i = 1; i < n; i++ )
	{
		double yi_1 = yi + h * f( xi, yi );
		yi = yi + h / 2.0 * ( f( xi, yi ) + f( xi + h, yi_1 ) );
		methods[ 1 ].push_back( yi );

		xi += h;
	}	
}

 
// Эйлер неявный
void Euler2()
{
	double yi = g_y0;
	double xi = 0.0f;

	methods[ 2 ].push_back( g_y0 );

	for( int i = 1; i < n; i++ )
	{
		double yi_1 = ( yi + h * f( xi, yi ) / 2.0 ) / ( 1.0 - h * f( xi + h, 1.0 ) / 2.0 );
		yi = yi + h / 2.0 * ( f( xi, yi ) + f( xi + h, yi_1 ) );
		methods[ 2 ].push_back( yi );

		xi += h;
	}	
}


//
void RungeKutta()
{
	double yi = g_y0;
	double xi = 0.0f;

	methods[ 3 ].push_back( g_y0 );

	for( int i = 1; i < n; i++ )
	{
		double k1 = h * f( xi, yi );
		double k2 = h * f( xi + h / 2.0, yi + k1 / 2.0 );
		double k3 = h * f( xi + h / 2.0, yi + k2 / 2.0 );
		double k4 = h * f( xi + h, yi + k3 );

		yi = yi + ( k1 + 2.0 * k2 + 2.0 * k3 + k4 ) / 6.0;

		methods[ 3 ].push_back( yi );

		xi += h;
	}		
}


//
void Koshi()
{
	double yi = g_y0;
	double xi = 0.0f;

	methods[ 4 ].push_back( g_y0 );

	for( int i = 1; i < n; i++ )
	{
		double yi_1 = yi + h / 2.0 * f( xi, yi );
		yi = yi + h * f( xi + h / 2.0, yi_1 );

		methods[ 4 ].push_back( yi );

		xi += h;
	}		
}


//
void Adams2Steps()
{
	double yi_1 = g_y0;
	double xi_1 = 0.0f;
	double yi = methods[ 3 ][ 1 ];
	double xi = xi_1 + h;

	methods[ 5 ].push_back( yi_1 );
	methods[ 5 ].push_back( yi );

	for( int i = 2; i < n; i++ )
	{
		double y = yi + h * ( 3.0 * f( xi, yi ) - f( xi_1, yi_1 ) ) / 2.0;

		yi_1 = yi;
		yi = y;

		methods[ 5 ].push_back( yi );

		xi_1 += h;
		xi += h;
	}		
}


//
void Adams3Steps()
{
	double yi_2 = g_y0;
	double xi_2 = 0.0f;
	double yi_1 = methods[ 3 ][ 1 ];
	double xi_1 = xi_2 + h;
	double yi = methods[ 3 ][ 2 ];
	double xi = xi_1 + h;

	methods[ 6 ].push_back( yi_2 );
	methods[ 6 ].push_back( yi_1 );
	methods[ 6 ].push_back( yi );

	for( int i = 3; i < n; i++ )
	{
		double y = yi + h * ( 23.0 * f( xi, yi ) - 16.0 * f( xi_1, yi_1 ) + 5.0 * f( xi_2, yi_2 ) ) / 12.0;

		yi_2 = yi_1;
		yi_1 = yi;
		yi = y;

		methods[ 6 ].push_back( yi );

		xi_2 += h;
		xi_1 += h;
		xi += h;
	}		
}


//
void Simpson()
{
	double yi = g_y0;
	double xi = 0.0f;

	methods[ 7 ].push_back( g_y0 );

	for( int i = 1; i < n; i++ )
	{
		double k1 = h * f( xi, yi );
		double k2 = h * f( xi + h / 2.0, yi + k1 / 2.0 );
		double k3 = h * f( xi + h, yi - k1 + 2 * k2 );

		yi = yi + k1 / 6.0f + k2 * 2.0f / 3.0f + k3 / 6.0f;

		methods[ 7 ].push_back( yi );

		xi += h;
	}		
}


// 0.1 e^(x (25.5-36.25 x+16.6667 x^2))
int main( int argc, char* argv[] )
{
	if( argc == 1 )
		return -1;

	FILE* file = fopen( "data", "w" );

	n = atoi( argv[ 1 ] );

	h = 1.0 / ( double )n;

	Euler();
	Euler2();
	Euler2W();
	RungeKutta();
	Koshi();
	Adams2Steps();
	Adams3Steps();
	Simpson();

	double xi = 0.0;

	for( size_t i = 0; i < n; i++ )
	{
		// точное решение
		double y = 0.1 * exp( xi * ( 25.5 - 36.25 * xi + 16.6667 * xi * xi ) );

		fprintf( file, "%0.6f %0.6f %0.6f %0.6f %0.6f %.06f %0.6f %0.6f %0.6f %0.6f\n", xi, methods[ 0 ][ i ], methods[ 1 ][ i ], methods[ 2 ][ i ], methods[ 3 ][ i ], methods[ 4 ][ i ], methods[ 5 ][ i ], methods[ 6 ][ i ], methods[ 7 ][ i ], y );

		xi += h;
	}

	fclose( file );

	return 0;
}
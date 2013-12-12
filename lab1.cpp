#include <stdio.h>
#include <math.h>
#include <functional>


// ф-ия из первой части лабы
double f1( double x )
{
	return sqrtf( 1.0f + 11.0f * cosf( x ) );
}


// из второй части лабы
double f2( double x )
{
	return 1.0f / ( 1.0f - x + x * x );
}


//
typedef std::function< double( double ) > Func;
typedef std::function< double( Func, int, double, double ) > Method;


// метод ср. прямоугольников
double AvgSquare( Func func, int n, double a, double b )
{
	double h = ( b - a ) / ( double )n;

	double sum = 0.0f;

	for( int i = 0; i < n; i++ )
	{
		double x = ( double )i * h;
		sum += func( x - h / 2.0f );
	}

	return h * sum;
}


// трапеций
double Trapeze( Func func, int n, double a, double b )
{
	double h = ( b - a ) / ( double )n;

	double sum = 0.0f;

	for( int i = 1; i < n; i++ )
	{
		double xi = a + ( double )i * h;

		sum += func( xi );
	}

	sum += func( a ) / 2.0f + func( b ) / 2.0f;

	return sum * h;
}



// симпосна
double Simpson( Func func, int n, double a, double b )
{
	double h = ( b - a ) / ( double )n;

	double sum = 0.0f;
	for( int i = 0; i < n; i += 2 )
	{
		double xi = ( double )i * h;
		sum += func( xi - h ) + 4.0f * func( xi ) + func( xi + h );
	}

	return sum * h / 3.0f;
}


// 
void ApplyMethod( Method method, int k )
{
	double h1 	= 		 method( f1, 10, 0.0f, 1.0f );
	double h2 	= 		 method( f1, 20, 0.0f, 1.0f );

	// правило Рунге
	double Rn = fabs( h2 - h1 ) / ( powf( 2.0f, ( double )k ) - 1.0f );

	printf( "%f %f\n", h1, Rn );
	printf( "%f\n", h2 );	
}


//
int main()
{
	// первая часть лабы
	ApplyMethod( AvgSquare, 2 );
	ApplyMethod( Trapeze, 2  );
	ApplyMethod( Simpson, 4 );

	// вторая
	FILE* file = fopen( "data", "w" );

	double f2Val = 2.0f * 3.1415926f / 3.0f / sqrt( 3.0f );
	for( int n = 2; n < 50; n += 2 )
	{
		double gAvgSq 		= fabs( AvgSquare( f2, n, 0.0f, 1.0f ) - f2Val);
		double gTrapeze 	= fabs( Trapeze( f2, n, 0.0f, 1.0f ) - f2Val );
		double gSimpson 	= fabs( Simpson( f2, n, 0.0f, 1.0f ) - f2Val );
		fprintf( file, "%d %.8f %.8f %.8f\n", n, gAvgSq, gTrapeze, gSimpson );
	}
	fclose( file );

	return 0;
}
#include <stdio.h>
#include <math.h>
#include <functional>


// ф-ия из первой части лабы
long double f1( long double x )
{
	return sqrtf( 1.0f + 11.0f * cosf( x ) );
}


// из второй части лабы
long double f2( long double x )
{
	return 1.0f / ( 1.0f - x + x * x );
}


//
typedef std::function< long double( long double ) > Func;
typedef std::function< long double( Func, int, long double, long double ) > Method;


// метод ср. прямоугольников
long double AvgSquare( Func func, int n, long double a, long double b )
{
	long double h = ( b - a ) / ( long double )n;

	long double sum = 0.0f;

	for( int i = 0; i < n; i++ )
	{
		long double x = ( long double )i * h;
		sum += func( x - h / 2.0f );
	}

	return h * sum;
}


// трапеций
long double Trapeze( Func func, int n, long double a, long double b )
{
	long double h = ( b - a ) / ( long double )n;

	long double sum = 0.0f;

	for( int i = 1; i < n; i++ )
	{
		long double xi = a + ( long double )i * h;

		sum += func( xi );
	}

	sum += func( a ) / 2.0f + func( b ) / 2.0f;

	return sum * h;
}



// симпосна
long double Simpson( Func func, int n, long double a, long double b )
{
	long double h = ( b - a ) / ( long double )n;

	long double sum = 0.0f;
	for( int i = 1; i < n; i += 2 )
	{
		long double xi = ( long double )i * h;
		sum += func( xi - h ) + 4.0f * func( xi ) + func( xi + h );
	}

	return sum * h / 3.0f;
}


// 
void ApplyMethod( Method method, int k )
{
	long double h1 	= 		 method( f1, 10, 0.0f, 1.0f );
	long double h2 	= 		 method( f1, 20, 0.0f, 1.0f );

	// правило Рунге
	long double Rn = fabs( h2 - h1 ) / ( powf( 2.0f, ( long double )k ) - 1.0f );

	printf( "%Lg %Lg\n", h1, Rn );
	printf( "%Lg\n", h2 );	
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

	long double f2Val = 2.0f * 3.1415926f / 3.0f / sqrt( 3.0f );
	for( int n = 2; n < 50; n += 2 )
	{
		long double gAvgSq 		= fabs( AvgSquare( f2, n, 0.0f, 1.0f ) - f2Val);
		long double gTrapeze 	= fabs( Trapeze( f2, n, 0.0f, 1.0f ) - f2Val );
		long double gSimpson 	= fabs( Simpson( f2, n, 0.0f, 1.0f ) - f2Val );
		fprintf( file, "%d %Lg %Lg %Lg\n", n, gAvgSq, gTrapeze, gSimpson );
	}
	fclose( file );

	return 0;
}

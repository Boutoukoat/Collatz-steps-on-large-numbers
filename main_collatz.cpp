// -----------------------------------------------------------------------
// Collatz step function calculator 
// -----------------------------------------------------------------------

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>

#include "collatz.h"
#include "bison.gmp_expr.h"

int main(int argc, char **argv)
{
	mpz_t n;
	mpz_init(n);

	for (int i = 1; i < argc; i++) {
		struct timespec ts1, ts2;
		unsigned long long steps;

		// Read an expression from the command line
		// Supported operators are +/-*^() with usual precedence.
                mpz_expression_parse(n, argv[i]);

		// gmp_printf("Test n=%Zd.\n",n);
		clock_gettime(CLOCK_REALTIME, &ts1);

		// Get the number of steps. Assuming the Collatz conjecture, it is not an infinite task. 
		steps = fastest_collatz(n);
		clock_gettime(CLOCK_REALTIME, &ts2);

		// Display the input number, the number of steps, and the time it took to run it, in milliseconds.
		double diff = ts2.tv_sec - ts1.tv_sec;
		diff *= 1e9;
		diff += ts2.tv_nsec;
		diff -= ts1.tv_nsec;
		diff /= 1e6;
		printf("f(%s)=%12llu, time=%12.3f msecs.\n", argv[i], steps, diff);
		fflush(stdout);
	}
	mpz_clear(n);

	return 0;
}


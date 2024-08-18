# Collatz-steps-on-large-numbers

Famous "Collatz 3n+1 iteration" can run on large numbers at a fraction of the textbook algorithms time.

From https://www.mersenneforum.org/showpost.php?p=634604&postcount=4, R.Gerbicz shows a serious speed-up 
changing the iteration complexity from O(log(n)^2) to O(log(n)^1.13). His GMP implementation can be 
seriously improved by a further O(magnitude factor), by moving out of GMP at the tail of the recursive
functions, and by a better control of temporary numbers and intermediate calculations.

To compile and run

sudo apt install libgmp3-dev
make

./collatz 10001
./collatz 10000 100000 1000000 10000000
./collatz 10000+1
./collatz "(5^3)*12"

As an example, this command 

./collatz 2^127-1 2^44497-1

displays

f(2^127-1)=        1660, time=       0.029 msecs.
f(2^44497-1)=      598067, time=      14.632 msecs.

There is an extremely small demo file "demo.py" to test with increasing numbers.

python demo.py
f(2^2-1)=           7, time=       0.001 msecs.
f(2^3-1)=          16, time=       0.000 msecs.
f(2^5-1)=         106, time=       0.000 msecs.
f(2^7-1)=          46, time=       0.000 msecs.
....
f(2^57885161-1)=   779044992, time=   29229.399 msecs.
f(2^74207281-1)=   998401306, time=   38161.140 msecs.
f(2^77232917-1)=  1039248803, time=   41093.431 msecs.
f(2^82589933-1)=  1111148968, time=   44362.609 msecs.
f(2^117385963-1)=  1579841291, time=   66985.520 msecs.
f(2^345678877-1)=  4651594256, time=  240021.755 msecs.

The largest number which can be tested depends on local memory, and on the gmp internal limits.


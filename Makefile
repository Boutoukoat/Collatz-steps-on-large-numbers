

GGG = g++ -O3 -march=native -fomit-frame-pointer -fexpensive-optimizations -flto

OBJ = main_collatz.o fastest_collatz.o bison.gmp_expr.o lex.gmp_expr.o


collatz: $(OBJ)
	$(GGG) -static -o collatz $(OBJ) -lgmp

main_collatz.o: main_collatz.cpp collatz.h
	$(GGG) -c -o main_collatz.o main_collatz.cpp

fastest_collatz.o: fastest_collatz.cpp collatz.h
	$(GGG) -c -o fastest_collatz.o fastest_collatz.cpp

bison.gmp_expr.o : bison.gmp_expr.tab.c bison.gmp_expr.h
	$(GGG) -c -o bison.gmp_expr.o bison.gmp_expr.tab.c

bison.gmp_expr.tab.c bison.gmp_expr.tab.h : parser.y
	bison -d parser.y

lex.gmp_expr.o : lex.gmp_expr.c
	$(GGG) -c -o lex.gmp_expr.o lex.gmp_expr.c

lex.gmp_expr.c : parser.l bison.gmp_expr.tab.h
	flex parser.l

clean:
	rm -f ./collatz $(OBJ)



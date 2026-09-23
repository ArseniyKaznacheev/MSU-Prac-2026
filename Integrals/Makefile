CC       = gcc
NASM     = nasm
CFLAGS   = -std=c99 -Wall -Wextra -O2 -m32
#ASMFLAGS = -f elf32
ASMFLAGS = -f win32 --prefix _
LDFLAGS  = -m32

.PHONY: all clean test

all: integral

funcs.o: funcs.asm
	$(NASM) $(ASMFLAGS) funcs.asm -o funcs.o

main.o: main.c
	$(CC) $(CFLAGS) -c main.c -o main.o

integral: main.o funcs.o
	$(CC) $(LDFLAGS) main.o funcs.o -o integral

#тестирование root и integral
test: integral
	@echo ________________________
	@echo -------test root------- 
	@echo ________________________
	./integral --test-root 1:2:0.0:1.5:0.0001:1.0
	./integral --test-root 5:6:0.0:2.0:0.0001:1.0
	./integral --test-root 3:4:-1.0:1.0:0.0001:0.0
	@echo ___________________________
	@echo -------test integral-------
	@echo ___________________________
	./integral --test-integral 1:0.0:1.0:0.0001:0.333333333
	./integral --test-integral 3:0.0:1.57079632:0.0001:1.0
	./integral --test-integral 5:0.0:1.0:0.0001:1.718281828

clean:
#	rm -f *.o integral
	del /Q *.o integral.exe 2>NUL
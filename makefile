simulator.so: simulator.c
	gcc -o $@ -fPIC -O3 -shared $^

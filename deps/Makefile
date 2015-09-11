libalnscore.so: alnscore.c
	$(CC) -dynamiclib -march=native -msse2 -O3 -o $@ $<
#	$(CC) -dynamiclib -msse2 -O0 -g -o $@ $<

clean:
	rm -rf libalnscore.so

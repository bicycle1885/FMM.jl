libalnscore.so: alnscore.c
	$(CC) -dynamiclib -march=native -msse2 -O3 -o $@ $<

.PHONY: clean
clean:
	rm -rf libalnscore.so

CC=gcc
CARGS=-Wall -O3
GMPLIB=-L/gmp_install/lib -lgmp
CMPHLIB=-L/usr/local/lib/libcmph.la -lcmph

karp-rabin:
	$(CC) $(CARGS) karp_rabin.c -o karp_rabin $(GMPLIB)

karp-rabin-clean:
	rm karp_rabin

hash-lookup:
	$(CC) $(CARGS) hash_lookup.c -o hash_lookup $(CMPHLIB)

hash-lookup-clean:
	rm hash_lookup

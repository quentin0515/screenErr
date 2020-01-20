all: screen htslib/lib/libhts.a

htslib/lib/libhts.a:
	cd htslib && autoheader && autoconf && ./configure --disable-s3 --disable-lzma --disable-bz2 --prefix=$(PWD)/htslib/ && make -j 4 && make install

OPTS=-O3
screen: screenErr.cpp htslib/lib/libhts.a
	g++ $(OPTS) -g -I htslib/ $^ -L htslib/lib -lhts -o $@  -Wl,-rpath,$(PWD)/htslib/lib -lz -lcurl -lpthread


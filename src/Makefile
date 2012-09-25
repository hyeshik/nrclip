CC=		gcc
SFMT_CFLAGS=	-O3 -finline-functions -fomit-frame-pointer -DNDEBUG \
		-fno-strict-aliasing --param max-inline-insns-single=1800 \
		-Wmissing-prototypes -Wall  -std=c99 -DSFMT_MEXP=19937
SFMT_SUBDIR=	SFMT-src-1.4
CLIPSIM_CFLAGS=	-O3 -Wall -DNDEBUG -I. -I${SFMT_SUBDIR}
CLIPSIM_LIBS=	-lm -lpthread -lz

.PHONY:	all

TARGETS=	crossfest SFMT.o

all: ${TARGETS}
	#time ./crossfest simdata-CLIP-2J3.error simdata-CLIP-2J3.reads 20 simdata-CLIP-2J3-
	python setup.py install --user
	PYTHONPATH=build/lib.linux-x86_64-2.7/ python test.py

crossfest: crossfest.c clipstatstree.c SFMT.o
	${CC} -o $@ ${CLIPSIM_CFLAGS} $+ ${CLIPSIM_LIBS}

SFMT.o: ${SFMT_SUBDIR}/SFMT.c
	${CC} ${SFMT_CFLAGS} -c -o $@ $<

clean:
	rm -f ${TARGETS}

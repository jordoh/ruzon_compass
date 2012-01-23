% separate targets for 32- and 64-bit platforms
% written originally for IRIX

COMPASSSRC = compassmex.c compass.c bs.c RGBLab.c emd.c
GREYSRC = greycompass.c RGBLab.c
GREYMATLABSRC = greycompassmex.c
GREYUNIXSRC = greycompassunix.c i_file.c i_main.c
BSSRC = bsmex.c RGBLab.c bs.c
CC = cc
CFLAGS =

all32:	compass32 bs32 grey32

all64:	compass64 bs64 grey64

compass32:	compassmex.mexsg

compassmex.mexsg:	$(COMPASSSRC) compass.h emd.h bs.h
	mex -O $(COMPASSSRC)

compass64:	compassmex.mexsg64

compassmex.mexsg64:	$(COMPASSSRC) compass.h emd.h bs.h
	mex -O $(COMPASSSRC)

grey32:	greycompassmex.mexsg greycompass

greycompassmex.mexsg:	$(GREYSRC) $(GREYMATLABSRC) greycompass.h
	mex -g $(GREYMATLABSRC) $(GREYSRC)

grey64:	greycompassmex.mexsg64 greycompass

greycompassmex.mexsg64:	$(GREYSRC) $(GREYMATLABSRC) greycompass.h
	mex -O $(GREYMATLABSRC) $(GREYSRC)

greycompass:	$(GREYSRC) $(GREYUNIXSRC) greycompass.h lib.h image.h
	cc -g -o greycompass $(GREYUNIXSRC) $(GREYSRC) -lm

bs32:	bsmex.mexsg

bsmex.mexsg:	$(BSSRC) bs.h 
	mex -O $(BSSRC)

bs64:	bsmex.mexsg64

bsmex.mexsg64:	$(BSSRC) bs.h
	mex -O $(BSSRC)


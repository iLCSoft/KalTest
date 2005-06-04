########################################################################	
##
## ---------------------------
##   Makefile for KalTest
## ---------------------------
##
## (Version Information)
##   Version: 0.2
##   Release: 1
##
## (Update Record)
##    2002/01/18  K.Hoshina     Derived from baby/src/Makefile
##    2002/10/21  K.Fujii       Cleanup.
##
## (Description)
##   In order to use this package you should first set some
##   environmental variables. Take a look at README in the
##   top directory.
##
## (Targets)
##	all       	: builds test programs
##	Makefiles 	: creates Makefiles.
##	clean     	: deletes *.o ...
##	distclean 	: deletes even test programs.
## 
########################################################################	

SHELL	= `which bash`
MFLAGS	=
CURRDIR	= .

SUBDIRS	= simple ct #cdc
SUBDIRS2 = hybrid

all:
	@case '${MFLAGS}' in *[ik]*) set +e;; esac; \
	for i in $(SUBDIRS); do \
	(cd $$i; echo ``making'' all ``in $(CURRDIR)/$$i...''; \
	xmkmf -a; \
	$(MAKE) $(MFLAGS)); \
	done
	@case '${MFLAGS}' in *[ik]*) set +e;; esac; \
	for i in $(SUBDIRS2); do \
	(cd $$i; echo ``making'' all ``in $(CURRDIR)/$$i...''; \
	$(MAKE) $(MFLAGS)); \
	done
	
Makefiles:
	@case '${MFLAGS}' in *[ik]*) set +e;; esac; \
	for i in $(SUBDIRS); do \
	(cd $$i; echo ``making'' Makefiles ``in $(CURRDIR)/$$i...''; \
	xmkmf); \
	done
	@case '${MFLAGS}' in *[ik]*) set +e;; esac; \
	for i in $(SUBDIRS2); do \
	(cd $$i; echo ``making'' Makefiles ``in $(CURRDIR)/$$i...''; \
	$(MAKE) $(MFLAGS) Makefiles); \
	done
	
clean: Makefiles
	@case '${MFLAGS}' in *[ik]*) set +e;; esac; \
	for i in $(SUBDIRS) $(SUBIRS2); do \
	(cd $$i; echo ``making'' clean ``in $(CURRDIR)/$$i...''; \
	$(MAKE) $(MFLAGS) clean); \
	done

distclean: clean
	@case '${MFLAGS}' in *[ik]*) set +e;; esac; \
	for i in $(SUBDIRS) $(SUBDIRS2); do \
	(cd $$i; echo ``making'' distclean ``in $(CURRDIR)/$$i...''; \
	$(MAKE) $(MFLAGS) distclean); \
	done;


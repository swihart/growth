LIBNAME=growth

OBJS=	carmasub.o cutil.o elliptic.o resid.o

$(LIBNAME)$(SHLIB_EXT):	$(OBJS)
	$(SHLIB_LD) $(SHLIB_LDFLAGS) -o $@ $(OBJS) $(FLIBS)

clean:
	@rm -f *.o *.$(SHLIB_EXT)

realclean: clean

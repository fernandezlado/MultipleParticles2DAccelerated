# 2D ACCELERATOR MAKEFILE

# --------------------------------------------------------------------------------
#	COMPILER			
# --------------------------------------------------------------------------------


FC=ifort


#DIRECTORIES


SRCDIR=src

OBJDIR=obj

MODDIR=mod

BINDIR=bin


# --------------------------------------------------------------------------------
# 	LIBRARIES
# --------------------------------------------------------------------------------


LIBS = -mkl


# --------------------------------------------------------------------------------
#	SOURCE FILES
# --------------------------------------------------------------------------------


SRC = \
	mkl_dfti.f90 \
	modMathConstants.f90 \
	modParameters.f90 \
	modIO.f90 \
	modFFT.f90 \
	modLinearAlgebra.f90 \
	modSpecialFunctions.f90 \
	modObstacle.f90 \
	modFarInteractions.f90 \
	modProjectionReferenceCell.f90 \
	modInterpolationReferenceCell.f90 \
	modCell.f90 \
	modForwardMap.f90 \
	main.f90


# --------------------------------------------------------------------------------
#	OBJECT FILES
# --------------------------------------------------------------------------------

# Replace .f90 extension by .o extension
OBJ_FILES = $(patsubst %.f90, %.o, $(SRC) )

# Add path to OBJDIR
OBJ = $(patsubst %, $(OBJDIR)/%, $(OBJ_FILES))

# --------------------------------------------------------------------------------
# 	Create Executable
# --------------------------------------------------------------------------------


#OUTPUT FILENAME

FILENAME = MultipleParticle

EXE = $(BINDIR)/$(FILENAME)


$(EXE): $(OBJ)
	$(FC) $(OBJ) $(LIBS) -o $(EXE)


$(OBJDIR)/%.o : $(SRCDIR)/%.f90
	$(FC) -c $< -o $@ -module $(MODDIR)/


clean:
	rm $(MODDIR)/*.mod
	rm $(OBJDIR)/*.o


cleanall:
	rm $(MODDIR)/*.mod
	rm $(OBJDIR)/*.o
	rm $(EXE)

run:
	./$(EXE)

debug: 
	$(FC) -ggdb $(patsubst %, $(SRCDIR)/%, $(SRC)) $(LIBS) -o $(EXE)_debug

ifeq (${FC},)
FC = gfortran
endif
ifeq (${FC},f77)
FC = gfortran
endif

FFLAGS =  -fPIC -cpp -fopenmp -O3 -ffree-line-length-none -std=f2008
DEBUGFLAGS = -g -O0 -fcheck=all -fbacktrace -cpp -fdec -fopenmp -ffree-line-length-none
#add -fpe0 to check for floating point errors (think lowLike also throws these harmlessly)
#MODOUT = -J$(OUTPUT_DIR)
LAPACKL = ${LAPACK_LINK}
CUBALIB = ../Cuba-1.4/libcuba.a

LDFLAGS = ${CUBALIB} ${LAPACKL}

all: lib


# Source-code directory
SRC_DIR = ../CosmoMC_BOSS/source

# Build directory
BUILD_DIR = build

BIN_DIR = bin
LIB_DIR = lib

# Debug build directory
DEBUG_BUILD_DIR = debug_build

# Objects
_OBJ = dummy_camb.o \
	   dummy_camb_constants.o \
	   MiscUtils.o \
	   MpiUtils.o \
	   StringUtils.o \
	   FileUtils.o \
	   IniObjects.o \
	   ObjectLists.o \
	   ObjectParamNames.o \
	   Matrix_utils_new.o \
	   settings.o \
	   IO.o \
	   GeneralTypes.o \
	   Interpolation.o \
	   likelihood.o \
	   CosmologyTypes.o \
	   CosmoTheory.o \
	   numrec.o \
	   quadpack.o \
	   rootfinding.o \
	   NLRSD_settings.o \
	   cubaint.o \
	   gRPT.o \
	   bias.o \
	   NLRSD_power.o \
	   Calculator_Cosmology.o \
	   fiducial.o \
	   Likelihood_Cosmology.o \
	   twopt_model.o \
	   twopt_data.o \
       twopt_likelihood.o \
	

# Add prefixes of build directory to objects
OBJ = $(addprefix $(BUILD_DIR)/,$(_OBJ))
DEBUG_OBJ = $(addprefix $(DEBUG_BUILD_DIR)/,$(_OBJ))

# ?
make_dirs = @mkdir -p $(@D)

all: wrapper

# Standard rules
lib: $(LIB_DIR)/libboss.a
wrapper: $(LIB_DIR)/libbias_wrapper.so
bin: $(BIN_DIR)/wrapper_driver

# Debugging rules
debug: FFLAGS += $(DEBUGFLAGS)
debug: $(LIB_DIR)/libboss.a

# Debugging rules
debug_bin: FFLAGS += $(DEBUGFLAGS)
debug_bin: $(BIN_DIR)/wrapper_driver_debug

# Rule to make object files
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.f90
	$(make_dirs)
	$(FC) -c -o $@ $< -J$(BUILD_DIR) $(LDFLAGS) $(FFLAGS)

# Rule to make debugging objects
$(DEBUG_BUILD_DIR)/%.o: $(SRC_DIR)/%.f90
	$(make_dirs)
	$(FC) -c -o $@ $< -J$(DEBUG_BUILD_DIR) $(LDFLAGS) $(FFLAGS)


# Rule to make HMx static library
$(LIB_DIR)/libboss.a: $(OBJ)
	@echo "\nBuilding static library.\n"
	$(make_dirs)
	$(AR) rc $@ $^

$(LIB_DIR)/libbias_wrapper.so: $(OBJ) bias_module.f90
	@echo "\nBuilding bias module wrapper.\n"
	@echo ${FC}
	$(make_dirs)
	$(FC) -shared -o $@ $^ -J$(BUILD_DIR) $(LDFLAGS) $(FFLAGS)

$(BIN_DIR)/wrapper_driver: $(OBJ) bias_module_test.f90
	@echo "\nBuilding executable.\n"
	$(make_dirs)
	$(FC) -o $@ $^ -J$(BUILD_DIR) $(LDFLAGS) $(FFLAGS)

$(BIN_DIR)/wrapper_driver_debug: $(DEBUG_OBJ) bias_module_test.f90
	@echo "\nBuilding executable.\n"
	$(make_dirs)
	$(FC) -o $@ $^ -J$(DEBUG_BUILD_DIR) $(LDFLAGS) $(FFLAGS)

# Rule to make cosmosis interface
#$(LIB_DIR)/libboss.so: $(SRC_DIR)/cosmosis_interface.f90
#	@echo "\nBuilding cosmosis interface.\n"
#	$(FC) $(FFLAGS) -shared -o $@ $^ -L$(LIB_DIR) -lboss $(LDFLAGS) -I$(BUILD_DIR) -J$(BUILD_DIR)

# Clean up
.PHONY: clean
clean:
	rm -f $(LIB_DIR)/libboss.a
	rm -f $(LIB_DIR)/libbias_wrapper.so
	rm -f $(BUILD_DIR)/*.o
	rm -f $(BUILD_DIR)/*.mod
	rm -f $(SRC_DIR)/*.mod
	rm -f $(DEBUG_BUILD_DIR)/*.o
	rm -f $(DEBUG_BUILD_DIR)/*.mod
	test -n "$(LIB_DIR)" && rm -rf $(LIB_DIR)/libbias_wrapper.so.dSYM/

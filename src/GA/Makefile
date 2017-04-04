
CXX=g++

EXE=elfin

SRC_DIR 		:= .
JUTIL_DIR 		:= ./jutil/src/
OBJ_DIR 		:= .obj
SRC_TREE 		:= $(shell find . -type d | grep -v './.obj\|\.dSYM\|^.$$')
OBJ_TREE 		:= $(addprefix $(OBJ_DIR)/, $(SRC_TREE))
BIN_DIR 		:= ./bin/
$(shell mkdir -p $(OBJ_TREE) $(BIN_DIR))

C_SRC 			:= $(shell find . -name '*.c')
CPP_SRC 		:= $(shell find . -name '*.cpp')

OBJS 			:= $(C_SRC:%.c=$(OBJ_DIR)/%.o) $(CPP_SRC:%.cpp=$(OBJ_DIR)/%.o)
DEPS 			:= $(C_SRC:%.c=$(OBJ_DIR)/%.d) $(CPP_SRC:%.cpp=$(OBJ_DIR)/%.d)

$(info Objects to be compiled: [${OBJS}])
#:

DEBUG=no
OMP=yes
TIMING=yes

ifeq ($(DEBUG), yes)
	DEBUG_FLAGS=-ggdb3
else
	DEBUG_FLAGS=
endif

ifeq ($(OMP), yes)
	ifeq ($(CXX), clang++)
		OMP_FLAGS=-openmp
	else
		OMP_FLAGS=-fopenmp -D_GLIBCXX_PARALLEL
	endif
else
	OMP_FLAGS=-D_NO_OMP
endif

ifeq ($(TIMING), yes)
	TIMING_FLAGS=-D_DO_TIMING
else
	TIMING_FLAGS=
endif

INCS 			:= -I./jutil/src/ -I.

ifeq ($(CXX), clang++)
	# For clang, these include directories vary from system to system
	# Find the one your gcc/g++ use
	INCS	 	+= -I/usr/local/Cellar/gcc/6.1.0/include/c++/6.1.0 \
				   -I/usr/local/Cellar/gcc/6.1.0/include/c++/6.1.0/x86_64-apple-darwin15.5.0 \
				   -I/usr/local/Cellar/gcc/6.1.0/include/c++/6.1.0/backward \
				   -I/usr/local/Cellar/gcc/6.1.0/lib/gcc/6/gcc/x86_64-apple-darwin15.5.0/6.1.0/include \
				   -I/usr/local/include \
				   -I/usr/local/Cellar/gcc/6.1.0/include \
				   -I/usr/local/Cellar/gcc/6.1.0/lib/gcc/6/gcc/x86_64-apple-darwin15.5.0/6.1.0/include-fixed
	# You also need libiomp (I got this from homebrew; your version may differ)
	INCS 		+= -I/usr/local/Cellar/libiomp/20150701/include/libiomp

	LD_FLAGS 	+= -Wl,-rpath,-L/usr/local/Cellar/gcc/6.1.0/lib/gcc/6 -lgomp

	ERR_FLAGS 	:=
	CPP_FLAGS 	+= -stdlib=libstdc++ 
else ifeq ($(CXX), g++)
	ERR_FLAGS	:= -fdiagnostics-color=always -fmax-errors=1
endif

OPT_FLAGS 		+= -O3
CPP_FLAGS 		+= -MMD -std=c++11 \
					$(OPT_FLAGS) $(DEBUG_FLAGS) $(OMP_FLAGS) $(TIMING_FLAGS) $(DEFS) $(INCS) $(EXTRA_FLAGS)

COMPILE 		:= $(CXX) $(CPP_FLAGS) $(ERR_FLAGS)


#
# start of rules
#

EXTS=c cpp
define make_rule
$(OBJ_DIR)/%.o: %.$1
	$$(COMPILE) -o $$@ -c $$< $(EXTRA_FLAGS)
endef
$(foreach EXT,$(EXTS),$(eval $(call make_rule,$(EXT))))

$(EXE): deleteTestObjs $(OBJS)
	$(COMPILE) $(OBJS) -o $(BIN_DIR)/$@ $(LD_FLAGS)

deleteTestObjs:
	rm -rf $(objToDelete)


test: $(EXE)
	./$(BIN_DIR)/$(EXE)

FORCE:
.PHONY: all clean

all: $(EXE)

clean: FORCE
	rm -rf $(BIN_DIR)/* $(OBJ_DIR)/* *.dSYM .DS_Store *.dec *.bin

-include $(DEPS)


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
	OMP_FLAGS=-fopenmp -D_GLIBCXX_PARALLEL
else
	OMP_FLAGS=-D_NO_OMP
endif

ifeq ($(TIMING), yes)
	TIMING_FLAGS=-D_DO_TIMING
else
	TIMING_FLAGS=
endif

INCS 			:= -I./jutil/src/ -I.

OPT_FLAGS 		:= -O3
CPPFLAGS 		:= -fdiagnostics-color=always -MMD -std=c++11 -fmax-errors=1\
					$(OPT_FLAGS) $(DEBUG_FLAGS) $(OMP_FLAGS) $(TIMING_FLAGS) $(DEFS) $(INCS)

COMPILE 		:= $(CXX) $(CPPFLAGS)


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
	$(COMPILE) $(LDFLAGS) $(OBJS) -o $(BIN_DIR)/$@ $(LDLIBS)

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

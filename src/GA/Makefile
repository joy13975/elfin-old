
CXX=g++

EXE=elfin

SRC_DIR 		:= .
JUTIL_DIR 		:= ./jutil/src/
OBJ_DIR 		:= .obj
SRC_TREE 		:= $(shell find . -type d | grep -v './.obj\|\.dSYM\|^.$$')
OBJ_TREE 		:= $(addprefix $(OBJ_DIR)/, $(SRC_TREE))
$(shell mkdir -p $(OBJ_TREE))

C_SRC 			:= $(shell find . -name '*.c')
CPP_SRC 		:= $(shell find . -name '*.cpp')

OBJS 			:= $(C_SRC:%.c=$(OBJ_DIR)/%.o) $(CPP_SRC:%.cpp=$(OBJ_DIR)/%.o)
DEPS 			:= $(C_SRC:%.c=$(OBJ_DIR)/%.d) $(CPP_SRC:%.cpp=$(OBJ_DIR)/%.d)

$(info Objects to be compiled: [${OBJS}])
#:

DEBUG=no
# OMP=yes

ifeq ($(DEBUG), yes)
	DEBUG_FLAG=-ggdb3
else
	DEBUG_FLAG=
endif

ifeq ($(GOMP), yes)
	GOMP_FLAG=-fopenmp
else
	GOMP_FLAG=
endif

INCS 			:= -I./jutil/src/

OPT_FLAGS 		:= -O3
CPPFLAGS 		:= -fdiagnostics-color=always -MMD -std=c++11 -fmax-errors=1\
					$(OPT_FLAGS) $(DEBUG_FLAG) $(GOMP_FLAG) $(DEFS) $(INCS)

COMPILE 		:= $(CXX) $(CPPFLAGS)


#
# start of rules
#

EXTS=c cpp
define make_rule
$(OBJ_DIR)/%.o: %.$1
	$$(COMPILE) -o $$@ -c $$<
endef
$(foreach EXT,$(EXTS),$(eval $(call make_rule,$(EXT))))

$(EXE): $(OBJS)
	$(COMPILE) $(LDFLAGS) $^ $(LDLIBS) -o $@

test: $(EXE)
	./$(EXE)

.PHONY: all clean

all: $(EXE)

clean:
	rm -rf $(EXE) $(OBJ_DIR)/* $(OBJ_DIR) *.dSYM .DS_Store *.dec *.bin

-include $(DEPS)
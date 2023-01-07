
#NOTES: 
# - $^ is the left side of the rule (the preconditions)
# - $@ is the right side of the rule (the name of the rule)
#
# - the function $(addprefix prefix, string_list) adds the prefix to all the 
#   elements in string_list
# - $(string_list:old_suffix=new_suffix) is very handy for going from .c files
#   to .o files (or similar).
# - $(patsubst src_regex,dest_regex, string_list) is a string substitution utility.
#   it substitute exploiting the 2 regexes all the elements of string_list
#

#STATIC_LIBS=./lib/libseal.a ./lib/libhdf5_hl_cpp.a ./lib/libhdf5_cpp.a ./lib/libhdf5_hl.a ./lib/libhdf5.a

#put source here
SRC_DIR=src

SRC = convolutionalLayer.cpp poolingLayer.cpp avgPoolingLayer.cpp fullyConnectedLayer.cpp squareLayer.cpp layer.cpp  sealencryptor.cpp utils.cpp cnnBuilder.cpp network.cpp H5Easy.cpp mainparams.cpp 
#your compiler
CXX = c++

#your flags
CXXFLAGS= -std=c++17 -I./$(SRC_DIR) -I./include -I./include/h5


LDFLAGS=-lseal -lpthread -Wl,-rpath,./lib -lz -ldl -lm
STATIC_LIBS=./lib/libseal.a -L./lib ./lib/libhdf5_cpp.a ./lib/libhdf5_hl_cpp.a  ./lib/libhdf5_hl.a ./lib/libhdf5.a

#executables will be saved here
BIN_DIR=bin
#object files will be saved here (for the reuse)
OBJ_DIR=obj

OBJECTS =$(SRC:.cpp=.o)
MAIN_NAME =main
MAIN_RULE = $(addprefix $(BIN_DIR)/, $(MAIN_NAME))
#all the object files excepts main object files
OBJECTS_PATH =$(addprefix $(OBJ_DIR)/,$(OBJECTS)) 
#all the object files
MAIN_OBJ_PATH=$(addprefix $(OBJ_DIR)/,$(MAIN_OBJ))

#if nothing is passed to the make command do this. It the starting point of the
# compilation
all: $(MAIN_RULE)

#build executables in bin_dir. Creates bin_dir and obj_dir if needed. 
# also compile objects files using prerequisites.
#   NB: it compiles also main object files that are not needed for the current
#       rule. A future fix will be to compile only the correct main object
$(MAIN_RULE): $(BIN_DIR) $(OBJ_DIR) $(OBJECTS_PATH) $(MAIN_OBJ_PATH) 
	$(CXX) $(OBJECTS_PATH) -o $@ $(CXXFLAGS) $(STATIC_LIBS) $(LDFLAGS)

#build objects taking source from src dir
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(CXX) $(CXXFLAGS) -c $^ -o $@ 

#create dirs
$(BIN_DIR):
	mkdir -p $(BIN_DIR)

$(OBJ_DIR):
	mkdir -p $(OBJ_DIR)

#remove all the object files
clean:
	rm -f $(OBJ_DIR)/*.o


#clean and build
full: clean all

CC      =	gcc
#CC      =   clang
CFLAGS  =	-Wall -O3 -Wno-unused-variable -Wno-unused-but-set-variable -Wno-unused-function
DFLAGS  =	-g -Wall
LIB     =	-lm -lz -lpthread
BIN_DIR =	./bin
SRC_DIR =   ./src

SOURCE  =	$(SRC_DIR)/abpoa.c $(SRC_DIR)/abpoa_align.c $(SRC_DIR)/abpoa_graph.c $(SRC_DIR)/simd_abpoa_align.c $(SRC_DIR)/simd_check.c $(SRC_DIR)/utils.c $(SRC_DIR)/abpoa_graph_visual.c
HEADER  =	$(SRC_DIR)/abpoa.h $(SRC_DIR)/abpoa_align.h $(SRC_DIR)/abpoa_graph.h $(SRC_DIR)/align.h $(SRC_DIR)/abpoa_graph_visual.h $(SRC_DIR)/kdq.h $(SRC_DIR)/kseq.h $(SRC_DIR)/ksort.h $(SRC_DIR)/simd_instruction.h $(SRC_DIR)/simd_abpoa_align.h $(SRC_DIR)/utils.h
OBJS    =	$(SRC_DIR)/abpoa_align.o $(SRC_DIR)/abpoa_graph.o $(SRC_DIR)/simd_abpoa_align.o $(SRC_DIR)/simd_check.o $(SRC_DIR)/utils.o $(SRC_DIR)/abpoa_graph_visual.o

# SIMD label
SSE41 			= __SSE4_1__
AVX2 			= __AVX2__
AVX512F 		= __AVX512F__
AVX512BW 		= __AVX512BW__

FLAG_SSE4       = -msse4
FLAG_AVX2       = -mavx2
FLAG_AVX512F    = -mavx512f
FLAG_AVX512BW   = -mavx512bw
SIMD_FLAG       =

.PHONY: all clean check
.SUFFIXES:.c .o

.c.o:
		$(CC) -c $(CFLAGS) $< -o $@

BIN     		= $(BIN_DIR)/abPOA
BPOALIB         = $(SRC_DIR)/libabpoa.a

SIMD_CHECK  	= $(BIN_DIR)/simd_check
GDB_DEBUG   	= $(BIN_DIR)/gdb_abPOA
SIMD_CHECK_D	= -D __CHECK_SIMD_MAIN__
DMARCRO 		= -D __DEBUG__

simd_flag := ${shell ./bin/simd_check 2> /dev/null}

ifeq ($(simd_flag), $(AVX512BW))
	SIMD_FLAG = $(FLAG_AVX512BW)
else ifeq ($(simd_flag), $(AVX512F))
	SIMD_FLAG = $(FLAG_AVX512F)
else ifeq ($(simd_flag), $(AVX2))
	SIMD_FLAG = $(FLAG_AVX2)
else ifeq ($(simd_flag), $(SSE41))
	SIMD_FLAG = $(FLAG_SSE4)
endif

all:		    $(BIN) 
abPOA:     		$(BIN)
gdb_abPOA: 		$(SOURCE) $(HEADER) $(GDB_DEBUG) 
libabpoa:        $(BPOALIB)


simd_check:$(SIMD_CHECK)
	$(shell ~/program/abPOA/bin/simd_check > /dev/null)

$(SIMD_CHECK):$(SRC_DIR)/simd_check.c $(SRC_DIR)/simd_instruction.h
	if [ ! -d $(BIN_DIR) ]; then mkdir $(BIN_DIR); fi
	$(CC) $(SIMD_CHECK_D) $(SRC_DIR)/simd_check.c -o $@

$(BIN):$(SRC_DIR)/abpoa.o $(BPOALIB)
	if [ ! -d $(BIN_DIR) ]; then mkdir $(BIN_DIR); fi
	$(CC) $(CFLAGS) $(SRC_DIR)/abpoa.o -o $@ -L$(SRC_DIR) -labpoa $(LIB) 

#TODO example.c

$(BPOALIB):$(OBJS)
	$(AR) -csr $@ $(OBJS)

$(SRC_DIR)/abpoa.o:$(SRC_DIR)/abpoa.c $(SRC_DIR)/abpoa.h $(SRC_DIR)/abpoa_graph.h $(SRC_DIR)/abpoa_align.h \
				  $(SRC_DIR)/abpoa_graph_visual.h $(SRC_DIR)/align.h $(SRC_DIR)/utils.h $(SRC_DIR)/simd_instruction.h
	$(CC) -c $(CFLAGS) $(SIMD_FLAG) $< -o $@

$(SRC_DIR)/simd_check.o:$(SRC_DIR)/simd_check.c $(SRC_DIR)/simd_instruction.h
	$(CC) -c $(CFLAGS) $(SIMD_FLAG) $< -o $@

$(SRC_DIR)/simd_abpoa_align.o:$(SRC_DIR)/simd_abpoa_align.c $(SRC_DIR)/abpoa_graph.h $(SRC_DIR)/abpoa_align.h $(SRC_DIR)/simd_instruction.h $(SRC_DIR)/utils.h
	$(CC) -c $(CFLAGS) $(SIMD_FLAG) $< -o $@

$(GDB_DEBUG): $(SOURCE) $(HEADER)
	if [ ! -d $(BIN_DIR) ]; then mkdir $(BIN_DIR); fi
	$(CC) $(DFLAGS) $(SIMD_FLAG) $(SOURCE) $(DMARCRO) -o $@ $(LIB)

clean:
	rm -f $(SRC_DIR)/*.[oa] $(BIN) $(SIMD_CHECK)

clean_debug:
	rm -f $(GDB_DEBUG)

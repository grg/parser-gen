-include Makefile.local

################################################################################
# Makefile Definitions
################################################################################


M_HOST := $(shell hostname)
PROJ_HOME := $(dir $(word $(words $(MAKEFILE_LIST)), $(MAKEFILE_LIST)))..
$(warning Warning: PROJ_HOME set to $(PROJ_HOME)) 


########### Generic Env Defs ############
#########################################


############# Design Defs ###############
#########################################

DESIGN_NAME:= parser
SRC_DIR := $(PROJ_HOME)/src
LIB_DIR := $(PROJ_HOME)/lib
BIN_DIR := $(PROJ_HOME)/bin
BUILD_DIR := $(PROJ_HOME)/build
INST_NAME := parser
MOD_NAME := parser
TOP_NAME := top
SYNTH_TOP_PATH := top.parser
VALID_DESIGN_NAME := yes
GENESIS_PROJECT_LIBS := $(LIB_DIR)/Perl5

ifndef SYN_CLK_PERIOD
  SYN_CLK_PERIOD := 2.5
  $(warning Warning: Using default clock period SYN_CLK_PERIOD=$(SYN_CLK_PERIOD))
endif


##################################################
# Parser parameters
##################################################

ifndef PARSE_GRAPH
  $(error Error: Parse graph not specified. Specify with PARSE_GRAPH=<parse_graph>)
endif

ifndef WORD_WIDTH
  WORD_WIDTH := 4
  $(warning Warning: Setting WORD_WIDTH=4)
endif

ifndef PROG_PARSER
  PROG_PARSER := 0
  $(warning Warning: Setting PROG_PARSER=0)
endif

ifndef EV_WIDTH
  ifeq ($(PROG_PARSER), 1)
    EV_WIDTH := 4096
  else
    EV_WIDTH := 4096
    #EV_WIDTH := auto
  endif
  $(warning Warning: Setting EV_WIDTH=$(EV_WIDTH))
endif

ifndef EV_INPUTS
  ifeq ($(PROG_PARSER), 1)
    EV_INPUTS := 12
  else
    EV_INPUTS := $(WORD_WIDTH)
  endif
  $(warning Warning: Setting EV_INPUTS=$(EV_INPUTS))
endif

ifeq ($(PROG_PARSER), 1)
  ifndef PROG_LOOKUP_WORDS
    PROG_LOOKUP_WORDS := 2
    $(warning Warning: Setting PROG_LOOKUP_WORDS=$(PROG_LOOKUP_WORDS))
  endif

  ifndef PROG_LOOKUP_WORD_WIDTH
    PROG_LOOKUP_WORD_WIDTH := 2
    $(warning Warning: Setting PROG_LOOKUP_WORD_WIDTH=$(PROG_LOOKUP_WORD_WIDTH))
  endif

  ifndef PROG_BUF_WORD_WIDTH
    # Set the default buffer word width to either 16 or the next power of two
    # >= the word width
    ifeq ($(shell expr $(WORD_WIDTH) \< 16),1)
      PROG_BUF_WORD_WIDTH := 16
    else
      PROG_BUF_WORD_WIDTH := $(shell perl -e 'use POSIX; print 2 ** ceil(log($(WORD_WIDTH)) / log(2));')
    endif
  endif

  ifndef MAX_RD_AMT
    MAX_RD_AMT := $(shell expr $(PROG_BUF_WORD_WIDTH) - 1)
    $(warning Warning: Setting MAX_RD_AMT=$(MAX_RD_AMT))
  endif
endif


PARSE_GRAPH_BASE := $(basename $(PARSE_GRAPH))
PARSE_GRAPH_JSON := $(PARSE_GRAPH_BASE).json 
GENESIS_PREREQS := $(PARSE_GRAPH_JSON)

PARSER_PARAMS := $(PARSER_PARAMS) ParseGraph=$(PARSE_GRAPH_JSON) WordWidth=$(WORD_WIDTH)
PARSER_PARAMS := $(PARSER_PARAMS) EnableProgParser=$(PROG_PARSER) ExtractVectorWidth=$(EV_WIDTH)

printnum = $(shell perl -e 'printf "%0$(2)d", $(1);')

ifeq ($(PROG_PARSER), 1)
  MAKE_TCAM_PARAMS := --lookups $(PROG_LOOKUP_WORDS) \
                      --lookup-width $(PROG_LOOKUP_WORD_WIDTH) \
                      --window $(PROG_BUF_WORD_WIDTH) \
                      --max-skip $(MAX_RD_AMT) \
                      --extract-bytes $(EV_INPUTS) \
		      --extract \
		      --no-first-lookup-at-zero \
		      --no-print-tcam \
		      --save-tcam
  PARSE_GRAPH_TCAM := $(PARSE_GRAPH_BASE)-w$(call printnum,$(PROG_BUF_WORD_WIDTH),3)
  PARSE_GRAPH_TCAM := $(PARSE_GRAPH_TCAM)-ms$(call printnum,$(MAX_RD_AMT),3)
  PARSE_GRAPH_TCAM := $(PARSE_GRAPH_TCAM)-l$(call printnum,$(PROG_LOOKUP_WORDS),2)
  PARSE_GRAPH_TCAM := $(PARSE_GRAPH_TCAM)-lw$(call printnum,$(PROG_LOOKUP_WORD_WIDTH),2)
  PARSE_GRAPH_TCAM := $(PARSE_GRAPH_TCAM)-flaz0-e1
  PARSE_GRAPH_TCAM := $(PARSE_GRAPH_TCAM)-eb$(call printnum,$(EV_INPUTS),2).tcam

  PARSER_PARAMS := $(PARSER_PARAMS) ProgLookupWords=$(PROG_LOOKUP_WORDS)
  PARSER_PARAMS := $(PARSER_PARAMS) ProgLookupWordWidth=$(PROG_LOOKUP_WORD_WIDTH)
  PARSER_PARAMS := $(PARSER_PARAMS) MaxRdAmt=$(MAX_RD_AMT)
  PARSER_PARAMS := $(PARSER_PARAMS) ProgBufWordWidth=$(PROG_BUF_WORD_WIDTH)
  PARSER_PARAMS := $(PARSER_PARAMS) NumExtractVectorInputs=$(EV_INPUTS)
  PARSER_PARAMS := $(PARSER_PARAMS) TCAMTableFile=$(PARSE_GRAPH_TCAM)

  GENESIS_PREREQS := $(GENESIS_PREREQS) $(PARSE_GRAPH_TCAM)
endif



############# For Genesis2 ##############
#########################################
# tile is the top of the pre-processed hierarchy
TOP_MODULE := $(TOP_NAME)


# list src folders and include folders for Genesis compiler
GENESIS_SRC := 	-srcpath $(SRC_DIR)/rtl			\
		-srcpath $(SRC_DIR)/rtl-old		\
		-srcpath $(SRC_DIR)/verif	

GENESIS_INC := 	-incpath $(SRC_DIR)/rtl			\
		-incpath $(SRC_DIR)/rtl-old		\
		-incpath $(SRC_DIR)/verif

# list src folders and include folders for Makefile to find (I hate duplicates!!!)
vpath %.svp $(GENESIS_SRC)
vpath %.vph $(GENESIS_INC)

# List source files
GENESIS_ENV :=		$(notdir $(wildcard $(SRC_DIR)/verif/*.svp))

GENESIS_DESIGN := 	$(notdir $(wildcard $(SRC_DIR)/rtl/*.svp)) \
			$(notdir $(wildcard $(SRC_DIR)/rtl-old/*.svp))

GENESIS_INPUTS :=	$(GENESIS_ENV) $(GENESIS_DESIGN) 

# List of file lists
GENESIS_INCLUDES := $(foreach dir,$(GENESIS_INC),$(wildcard $(dir)/*.vph))

# debug level
GENESIS_DBG_LEVEL := 0

# List of generated verilog files
GENESIS_VLOG_LIST := genesis_vlog.vf
GENESIS_VLOG_SYNTH_LIST := genesis_vlog.synth.vf
GENESIS_VLOG_VERIF_LIST := genesis_vlog.verif.vf

# Input xml program
#ifndef GENESIS_CFG_XML
#  GENESIS_CFG_XML := 	$(SRC_DIR)/configs/empty_cfg.xml
#endif
#$(warning WARNING: GENESIS_CFG_XML set to $(GENESIS_CFG_XML))

# xml hierarchy file
ifndef GENESIS_HIERARCHY
  GENESIS_HIERARCHY := $(MOD_NAME).xml
endif
$(warning WARNING: GENESIS_HIERARCHY set to $(GENESIS_HIERARCHY))

# For more Genesis parsing options, type 'Genesis2.pl -help'
#        [-parse]                    ---   should we parse input file to generate perl modules?
#        [-sources|srcpath dir]      ---   Where to find source files
#        [-includes|incpath dir]     ---   Where to find included files
#        [-input file1 .. filen]     ---   List of top level files
#                                    ---   The default is STDIN, but some functions
#                                    ---   (such as "for" or "while")
#                                    ---   may not work properly.
#        [-perl_modules modulename]  ---   Additional perl modules to load
GENESIS_PARSE_FLAGS := 	-parse $(GENESIS_SRC) $(GENESIS_INC)

ifneq ($(strip $(GENESIS_INPUTS)),)
  GENESIS_PARSE_FLAGS := $(GENESIS_PARSE_FLAGS) -input $(GENESIS_INPUTS)
endif

ifneq ($(strip $(PARSER_PARAMS)),)
  space :=
  space +=
  GENESIS_PARAMETERS = $(subst  $(space), --parameter top.parser., $(strip $(PARSER_PARAMS)))
endif

ifneq ($(strip $(VERIF_PARAMS)),)
  space :=
  space +=
  GENESIS_PARAMETERS = $(subst  $(space), --parameter top., $(strip $(VERIF_PARAMS)))
endif
#GENSESIS_SYNTHESIS_OVERRIDE := -parameter top.GenMode=Synth

# For more Genesis parsing options, type 'Genesis2.pl -help'
#        [-generate]                 ---   should we generate a verilog hierarchy?
#        [-top topmodule]            ---   Name of top module to start generation from
#        [-depend filename]          ---   Should Genesis2 generate a dependency file list? (list of input files)
#        [-product filename]         ---   Should Genesis2 generate a product file list? (list of output files)
#        [-hierarchy filename]       ---   Should Genesis2 generate a hierarchy representation tree?
#        [-xml filename]             ---   Input XML representation of definitions
GENESIS_GEN_FLAGS :=	-gen -top $(TOP_MODULE)					\
			-depend depend.list					\
			-product $(GENESIS_VLOG_LIST)				\
			-hierarchy $(GENESIS_HIERARCHY)                		\
			-debug $(GENESIS_DBG_LEVEL)				\
			-synthtop $(SYNTH_TOP_PATH)                             
ifdef GENESIS_CFG_XML
  GENESIS_GEN_FLAGS := $(GENESIS_GEN_FLAGS)					\
			-xml $(GENESIS_CFG_XML)					
endif

GENESIS_GEN_FLAGS := $(GENESIS_GEN_FLAGS)					\
			$(GENESIS_PARAMETERS)




############### For Verilog (VCS) ##########
############################################

# Prevent simulating with clock frequencies that are too small
MIN_SIM_CLK_PERIOD := 0.02

SIM_CLK_PERIOD := $(SYN_CLK_PERIOD)
ifeq ($(origin SIM_CLK_PERIOD),file)
  ifeq ($(shell expr $(SIM_CLK_PERIOD) '<' $(MIN_SIM_CLK_PERIOD)),1)
    SIM_CLK_PERIOD := $(MIN_SIM_CLK_PERIOD)
  endif
endif

COMPILER := vcs

EXECUTABLE := ./simv

VERILOG_ENV :=		 

VERILOG_DESIGN ?=	

VERILOG_FILES :=  	$(VERILOG_ENV)	$(VERILOG_DESIGN)					

VERILOG_LIBS := 	

ifndef SYNOPSYS
  $(error ERROR: SYNOPSYS environment variable is undefined )
endif

ifeq ($(shell uname -p), x86_64)
    VERILOG_CMP_FULL64 ?= -full64
endif


# "-sverilog" enables system verilog
# "+lint=PCWM" enables linting error messages
# "+libext+.v" specifies that library files (imported by the "-y" directive) ends with ".v"
# "-notice" used to get details when ports are coerced to inout
# "-full64" for 64 bit compilation and simulation
# "+v2k" for verilog 2001 constructs such as generate
# "-timescale=1ns/1ns" sets the time unit and time precision for the entire design
# "+noportcoerce" compile-time option to shut off the port coercion for the entire design
# "-top topModuleName" specifies the top module
# "-f verilogFiles.list" specifies a file that contains list of verilog files to compile
# "+vcs+saif_libcell" 
# "-y $SYNOPSYS/dw/sim_ver" include for designware macros


VERILOG_COMPILE_BASE_FLAGS := 	-sverilog 					\
				+vcs+lic+wait					\
				+cli 						\
				+lint=PCWM					\
				+lint=UV					\
				+lint=IRIMW					\
				+lint=FIAAS					\
				+lint=CAWM-L					\
				+lint=ERASM-L					\
				-ova_lint -ova_lint_magellan			\
				+warn=all					\
				-Mupdate					\
				-picarchive -Marchive				\
				+libext+.v					\
				+libext+.inc					\
				-notice						\
				$(VERILOG_CMP_FULL64)  				\
				+v2k						\
				-debug_pp					\
				-timescale=1ps/1ps				\
				+noportcoerce         				\
				-ld $(J_CC)          				\
				-top $(TOP_MODULE)				\
				-parallel+show_features				\
				-parallel+sva=4                                 \
				-parallel+saif=1                                \
				+vcs+saif_libcell                               \
				$(SYN_LIB_VERILOG_PAD_FILES) $(VERILOG_FILES) $(VERILOG_LIBS) $(VERILOG_DW_MISC_FILES)

VERILOG_COMPILE_FLAGS := 	$(VERILOG_COMPILE_BASE_FLAGS)                   \
				-f $(GENESIS_VLOG_LIST)

VERILOG_COMPILE2_FLAGS := 	$(VERILOG_COMPILE_BASE_FLAGS)                   \
				-f $(GENESIS_VLOG_VERIF_LIST)


VCS_SAIF_FILE = $(TOP_NAME).saif 

# +VCS_SAIF_FILE=top.saif
# "+vpdbufsize+100" limit the internal buffer to 100MB (forces flushing to disk)
# "+vpdports" Record information about ports (signal/in/out)
# "+vpdfileswitchsize+1000" limits the wave file to 1G (then switch to next file)
ORIG_VERILOG_SIMULATION_FLAGS := $(VERILOG_SIMULATION_FLAGS)
#VERILOG_SIMULATION_FLAGS := 	$(VERILOG_SIMULATION_FLAGS) 			\

set_verilog_simulation_flags := VERILOG_SIMULATION_FLAGS := \
				$$(ORIG_VERILOG_SIMULATION_FLAGS) 		\
				+vcs+lic+wait					\
				+vpdbufsize+100					\
				+vpdfileswitchsize+100				\
				+vcs+saif_libcell                               \
				+VCS_SAIF_FILE=$$(VCS_SAIF_FILE)                 \
                                +clk_period=$$(SIM_CLK_PERIOD)
$(eval $(set_verilog_simulation_flags))


################################################################################
################ Makefile Rules
################################################################################
#default rule: 
all: $(EXECUTABLE)

# Header processing rules:
#############################
$(PARSE_GRAPH_JSON): $(PARSE_GRAPH)
	@echo Building JSON: $<
	$(BIN_DIR)/make_json.py $< $@

$(PARSE_GRAPH_TCAM): $(PARSE_GRAPH)
	@echo Building TCAM configuration: $<
	$(BIN_DIR)/make_tcam.py $(MAKE_TCAM_PARAMS) $<


# Genesis2 rules:
#####################
# Genesis2 Generate:
# This is the rule to activate Genesis2 generator to generate verilog 
# files (_unqN.v) from the generator program (.svp files).
# Use "make PARSE=<genesis2_parse_flags> GEN=<genesis2_gen_flags>" to add elaboration time flags
LOG_GEN := 2>&1 | tee --append gen2.log ; perl $(PROJ_HOME)/bin/checkRun.pl gen2.log
ifneq ($(GENESIS_PROJECT_LIBS),)
  GENESIS_ENV = GENESIS_PROJECT_LIBS=$(GENESIS_PROJECT_LIBS)
endif
$(GENESIS_VLOG_SYNTH_LIST): $(GENESIS_VLOG_LIST)
$(GENESIS_VLOG_VERIF_LIST): $(GENESIS_VLOG_LIST)
$(GENESIS_VLOG_LIST): $(GENESIS_INPUTS) $(GENESIS_CFG_XML) $(GENESIS_INCLUDES) $(GENESIS_PREREQS)
	@echo ""
	@echo Making $@ because of $?
	@echo ==================================================
	@hostname > hostname.log
	$(GENESIS_ENV) Genesis2.pl $(GENESIS_PARSE_FLAGS) $(PARSE) $(GENESIS_GEN_FLAGS) $(GEN) $(LOG_GEN)


# phony rules for partial compilation process
.PHONY: genesis_clean gen gen_syn

gen: $(GENESIS_VLOG_LIST)

LOG_GEN_SYN := 2>&1 | tee --append gen2_syn.log ; perl $(PROJ_HOME)/bin/checkRun.pl gen2_syn.log
gen_syn: $(GENESIS_INPUTS) $(GENESIS_CFG_XML)
	@echo ""
	@echo Making $@ because of $?
	@echo ==================================================
	rm -f *.v 
	$(GENESIS_ENV) Genesis2.pl $(GENESIS_PARSE_FLAGS) $(PARSE) $(GENESIS_GEN_FLAGS) $(GENSESIS_SYNTHESIS_OVERRIDE) $(GEN) $(LOG_GEN_SYN)

genesis_clean: 
	@echo ""
	@echo Cleanning previous runs of Genesis
	@echo ===================================
	\rm -rf $(GENESIS_VLOG_LIST) $(GENESIS_VLOG_SYNTH_LIST) $(GENESIS_VLOG_VERIF_LIST)
	if test -f "genesis_clean.cmd"; then    \
		$(SHELL) genesis_clean.cmd;       \
	fi




# VCS rules:
#####################
# compile rules:
# use "make COMP=+define+<compile_time_flag[=value]>" to add compile time flags
.PHONY: comp
comp: $(EXECUTABLE)

COMP_RETRIES := 5

LOG_COMP :=  2>&1 | tee --append comp.log
$(EXECUTABLE):	$(VERILOG_FILES) $(GENESIS_VLOG_LIST) 
	@echo ""
	@echo Making $@ because of $?
	@echo ==================================================
	@DONE=0 ; COUNT=0 ; \
	until [ "$$DONE" -eq 1 ] ; do \
		COUNT=`expr $$COUNT + 1` ; \
		if [ "$$COUNT" -eq $(COMP_RETRIES) ] ; then \
			DONE=1 ; \
		fi ; \
		$(COMPILER)  $(VERILOG_COMPILE_FLAGS) $(COMP) $(LOG_COMP) ; \
		if grep -q "Segmentation fault" comp.log ; then \
			mv comp.log comp.$${COUNT}.log ; \
		else \
			DONE=1 ; \
		fi ; \
	done
	@grep -v "Segmentation fault" comp.log && perl $(PROJ_HOME)/bin/checkRun.pl comp.log 

.PHONY: comp2

LOG_COMP_SYN :=  2>&1 | tee --append comp_syn.log
comp2: 
	@echo ""
	@echo Making $@ because of $?
	@echo ==================================================
	@DONE=0 ; COUNT=0 ; \
	until [ "$$DONE" -eq 1 ] ; do \
		COUNT=`expr $$COUNT + 1` ; \
		if [ "$$COUNT" -eq $(COMP_RETRIES) ] ; then \
			DONE=1 ; \
		fi ; \
		$(COMPILER)  $(VERILOG_COMPILE2_FLAGS) $(COMP) phys_smpl/results/$(MOD_NAME)_unq1.mapped.v $(SYN_LIB_VERILOG_FILES) $(LOG_COMP_SYN) ; \
		if grep -q "Segmentation fault" comp_syn.log ; then \
			mv comp_syn.log comp_syn.$${COUNT}.log ; \
		else \
			DONE=1 ; \
		fi ; \
	done
	@grep -v "Segmentation fault" comp_syn.log && perl $(PROJ_HOME)/bin/checkRun.pl comp_syn.log 


# Simulation rules:
#####################
# use "make run RUN=+<runtime_flag[=value]>" to add runtime flags
.PHONY: run run2 run_wave

LOG_RUN := ; perl $(PROJ_HOME)/bin/checkRun.pl simv.log 
LOG_RUN2 := ; perl $(PROJ_HOME)/bin/checkRun.pl simv2.log 
SIM_LOG := -l $(EXECUTABLE).log
SIM_LOG2 := -l $(EXECUTABLE)2.log
run_wave: $(EXECUTABLE)
	@echo ""
	@echo Now Running $(EXECUTABLE) with wave
	@echo ==================================================
	$(EXECUTABLE) +wave $(VERILOG_SIMULATION_FLAGS) $(SIM_LOG) $(RUN) $(LOG_RUN)

run: $(EXECUTABLE)
	@echo ""
	@echo Now Running $(EXECUTABLE)
	@echo ==================================================
	@echo with clock period $(SIM_CLK_PERIOD)
	$(EXECUTABLE) $(VERILOG_SIMULATION_FLAGS) $(SIM_LOG) $(RUN) $(LOG_RUN)

run-gui: $(EXECUTABLE)
	@echo ""
	@echo Now Running $(EXECUTABLE)
	@echo ==================================================
	@echo with clock period $(SIM_CLK_PERIOD)
	$(EXECUTABLE) -gui $(VERILOG_SIMULATION_FLAGS) $(SIM_LOG) $(RUN) $(LOG_RUN)

run2: $(EXECUTABLE)
	@echo ""
	@echo Now Running $(EXECUTABLE)
	@echo ==================================================
	@echo with clock period $(SIM_CLK_PERIOD)
	$(EXECUTABLE) $(VERILOG_SIMULATION_FLAGS) $(SIM_LOG2) $(RUN) $(LOG_RUN2)


# SIM
########################
.PHONY: sim
sim: gen comp run



# Cleanup rules:
#####################
.PHONY: clean cleanall 
clean:  genesis_clean
	@echo ""
	@echo Cleanning old files, objects, logs and garbage
	@echo ==================================================
	\rm -rf $(EXECUTABLE) 
	\rm -rf csrc
	\rm -rf *.daidir
	\rm -rf *.log
	\rm -rf *.pvl
	\rm -rf *.syn
	\rm -rf *.mr
	\rm -rf *.pvk
	\rm -rf *.flc
	\rm -rf ucli.key
	\rm -rf *~
	\rm -rf $(TOP_NAME).saif
	\rm -rf $(EXECUTABLE) $(EXECUTABLE).vdb


cleanall: clean clean_dc clean_icc clean_reports
	\rm -rf DVE*
	\rm -f  *_unq1.v_BACKUP
	\rm -rf vcdplus.vpd
	\rm -f *.v
	\rm -f *.pm

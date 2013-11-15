Parser generator
================

This project generates packet parsers for use in network devices such as
switches and routers. It generates bot fixed and programmable parsers. Fixed
parsers use a parse graph that is chosen at generation time, while
programmable parsers use a parse graph that is chosen at run time.

The generator was originally created to fascilitate exploration of the parser
design space. More information can be found in *Design Principles for Packet
Parsers* by Glen Gibb et al. (See references.)

The generator is *not* the same version used to produce results for the paper.
This version offers fewer configurable parameters in order to make the code
easier to understand and modify.


Directory structure
===================

    |
    +-- bin         - scripts used during generation
    |
    +-- build       - target directory for generated files
    |
    +-- examples    - example parse graphs
    |
    +-- lib         - libraries used by the various scripts
    |   |
    |   +-- Perl5   - Perl libraries
    |
    +-- setup       - environment configuration
    |
    +-- src         - Genesis source files
        |
        +-- rtl     - RTL for the parser
        |
        +-- verif   - testbench and associated files


Prerequisites
=============
The parser generator uses the Genesis chip generator to perform generation.
Information about obtaining Genesis can be found here:
http://genesis2.stanford.edu/

The generator produces synthesizable SystemVerilog. Simulation, synthesis, and
place-and-route require the use of third party tools. The tools used for
testing during development were:
 * Synopsys VCS Simulator
 * Synopsys Design Compiler
 * Synopsys IC Compiler
The generated code should work with any tool that supports SystemVerilog.


Usage
=====

A parse graph is required to use the parser generator. The logic in a fixed
parser is generated specifically for the chosen parse graph. The logic in a
programmable parser is independent of the parse graph; however, the parse
graph is used to build the test bench and to populate the parse state table
(emulated TCAM/RAM) module.


Parse graphs
------------

A parse graph is a directed acyclic graph that specifies the types and
permitted orderings of headers within a network. The parser identifies a path
within the parse graph for each received packet.

Parse graphs are described via text files that list a set of properties for
each header. Example parse graphs are found in the `examples` directory.
See README.parse-graphs.md for a description of the format.


Environment setup
-----------------

The Makefile must be able to execute Genesis and the simulator/synthesis
tools. The script in the setup directory provides an example of a bash
configuration to set paths appropriately. Customize this to your needs.
`build/Makefile` must be customized for the simulator and synthesis tools you
use.


Generating a parser
-------------------

Parsers are generated in the build directory. To build a parser, do:

    cd build
    make <TARGETS> <PARSER_PARAMS> <VERIF_PARAMS>

Valid targets are:
 * `clean` - clean the build directory of generated files.
 * `gen` - generate a new parser. Synthesizable SystemVerilog is placed in
   `genesis_synth` and the test code is place in `genesis_verif`.
 * `comp` - compile the generated design with the simulator.
 * `run` - run the compiled simulation.
 * `sim` - combination of `gen`, `run`, and `sim`.

`PARSER_PARAMS` allows parser parameters to be specified to the generator.
Parser parameters are specified as:

    param1=val1 param2=val2 .. paramN=valN

Parser parameters are:
 * `PARSE_GRAPH`: the parse graph file to use. __(Required.)__
 * `WORD_WIDTH`: input width of parser (bytes).
 * `PROG_PARSER`: generate a programmable parser? (0 = fixed, 1 =
   programmable.)
 * `EV_WIDTH`: width of the packet header vector (bits).

Additional parameters for programmable parsers:
 * `PROG_LOOKUP_WORDS`: number of inputs to the parser state table.
 * `PROG_LOOKUP_WORD_WIDTH`: width of each input to the state table.
 * `PROG_BUF_WORD_WIDTH`: width of internal buffer.
 * `EV_INPUTS`: number of extract vector inputs.
 * `MAX_RD_AMT`: maximum number of bytes to advance in a single cycle.


`VERIF_PARAMS` determines how the test bench/verification code is generated.
Verification parameters are specified as:

    VERIF_PARAMS="param1=val1 param2=val2 .. paramN=valN"

Verification parameters are:
 * `MaxPkts`: the maximum number of packets for the testbench to generate.
 * `RandomPktData`: fill the test packets with random data for fields not used
   to determine the next header type. (1 = use random data, 0 = no random
   data).
 * `SeqPktData`: fill the test packets with sequential data for fields not
   used to determine the next header type. (1 = use sequential data, 0 = do
   not use sequential data).

_`RandomPktData` and `SeqPktData` cannot both be used at the same time. The
sequential option overrides the random option if both are specified._

**Note: Generating a parser will overwrite any previously generated parsers.**

### Usage examples:

Generate a fixed parser for the parse graph `examples/graph-simple.txt` using
a word width of 4B:

    make gen PARSER_PARAMS="ParseGraph=../examples/graph-simple.txt \
                            WordWidth=4"

Generate a programmable parser for the parse graph
`examples/graph-enterprise.txt` using a word width of 8B:

    make gen PARSER_PARAMS="ParseGraph=../examples/graph-enterprise.txt \
                            WordWidth=8 \
                            EnableProgParser=1"


Verification
============

The test bench verifies the parser by sending in a sequence of test packets to
test each path in the parse graph (or a subset of paths determined by the
`MaxPkts` verification parameter). The test bench verifies that the correct
packet header vector is output for each input file.


References
==========

* _Design Principles for Packet Parsers._ Glen Gibb, George Varghese, Mark
Horowitz, and Nick McKeown:
http://dl.acm.org/citation.cfm?id=2537857.2537860&coll=DL&dl=GUIDE&CFID=222668031&CFTOKEN=52789323
* _Reconfigurable Hardware for Software-Defined Networks, Chapter 2._ Glen Gibb:
http://purl.stanford.edu/ns046rz4288

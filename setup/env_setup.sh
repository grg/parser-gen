# Bourne shell script to initialize the environment for the packet generator
#
# Note: this script must be loaded using source. e.g.: "source env_setup.sh"

# Variable that points at the top-level parser generator directory
PARSER_GEN=

if [ "$PARSER_GEN" == "" ] ; then
	echo $BASH_SOURCE: Error: please edit the script and set PARSER_GEN
fi

# Set the project Perl path
export GENESIS_PROJECT_LIBS=$PARSER_GEN/lib/Perl5
export PYTHONPATH=$PYTHONPATH:$PARSER_GEN/lib/python

# Set up other tools here, such as Genesis, simulator, synthesis, and backend

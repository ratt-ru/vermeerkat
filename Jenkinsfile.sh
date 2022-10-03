set -e

echo "----------------------------------------------"
echo "$JOB_NAME build $BUILD_NUMBER"
WORKSPACE_ROOT="$WORKSPACE/$BUILD_NUMBER"
echo "Setting up build in $WORKSPACE_ROOT"
TEST_OUTPUT_DIR_REL=testcase_output
TEST_OUTPUT_DIR="$WORKSPACE_ROOT/$TEST_OUTPUT_DIR_REL"
TEST_DATA_DIR="$WORKSPACE/../../../test-data"
export TEST_DATA_DIR
export TEST_OUTPUT_DIR
PROJECTS_DIR_REL="projects"
PROJECTS_DIR=$WORKSPACE_ROOT/$PROJECTS_DIR_REL
mkdir $TEST_OUTPUT_DIR
echo "----------------------------------------------"
echo "\nEnvironment:"
df -h .
echo "----------------------------------------------"
cat /proc/meminfo
echo "----------------------------------------------"

export SINGULARITY_PULLFOLDER=$PROJECTS_DIR/.stimela_images 
mkdir $PROJECTS_DIR/.stimela_images

cd $TEST_OUTPUT_DIR
virtualenv $PROJECTS_DIR/venv -p python3.6
source $PROJECTS_DIR/venv/bin/activate
python3.6 -m pip install -U pip setuptools wheel
python3.6 -m pip install $PROJECTS_DIR/VermeerKAT
python3.6 -m pip install nose

stimela pull -d
function finish {
    rm -rf $SINGULARITY_PULLFOLDER
}
trap finish EXIT

python3.6 -m nose $PROJECTS_DIR/VermeerKAT/tests/acceptance_test.py

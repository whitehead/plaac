#!/bin/bash

SCRIPT_DIR=`dirname $0`
cd $SCRIPT_DIR

# install Python's setuptools-scm to get version information
VENV=.build-venv

if [ ! -d "$VENV" ]; then
    echo "Creating Python build environment..."
    python3 -m venv "$VENV"
fi

. "$VENV/bin/activate"

python -m pip install --upgrade pip >/dev/null
python -m pip install -q "setuptools-scm>=8,<9"

# now get the package version

echo "getting version"

PLAAC_VERSION=$(python - <<'EOF'
from setuptools_scm import get_version
print(get_version(
    root="..",
    version_scheme="guess-next-dev",
    local_scheme="node-and-date",
    fallback_version="1.1.0.dev0",
))
EOF
)

echo "ok (Version: ${PLAAC_VERSION})"

echo cleanup
rm -rf target
mkdir target

# generate version using the tag from github 
mkdir -p target/generated
cat > target/generated/Version.java <<EOF
public class Version {
    public static final String VERSION = "$PLAAC_VERSION";
}
EOF

echo build
javac --release 11 -d target src/*.java target/generated/Version.java

echo package
cp src/mainClass target/.
cp -r src/util target/.
cd target
jar cmf mainClass plaac.jar *
chmod +x plaac.jar

echo "ok (see target/plaac.jar)"

# also build PLAAC docs
echo "build-docs"

cd -
# generate documentation for command-line parameters, and dotfile for default hmm
java -jar target/plaac.jar -s -d -i dummy -h target/hmm_default.txt > target/plaac_headers.txt
python3 ./build_docs.py target/plaac_headers.txt > target/_plaac_headers.haml

echo "ok (see target/_plaac_headers.haml)"

if ! type "dot" > /dev/null; then
  echo "dot not found, skipping genereration of hmm_default.png" 
else
  dot target/hmm_default.txt -Tpng > target/hmm_default.png
  echo "ok (see target/hmm_default.png)" 
fi

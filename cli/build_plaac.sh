#!/bin/bash

SCRIPT_DIR=`dirname $0`
cd $SCRIPT_DIR

# now get the package version

if [ -n "${PLAAC_VERSION:-}" ]; then
    echo "using supplied version: ${PLAAC_VERSION}"
else
    echo "getting version"
    PLAAC_VERSION=$(./get_plaac_version.sh)
fi

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

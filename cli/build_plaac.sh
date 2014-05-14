
SCRIPT_DIR=`dirname $0`
cd $SCRIPT_DIR

echo cleanup
rm -rf target
mkdir target

echo build
javac -source 1.5 -target 1.5 -d target src/*.java

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
./build_docs.py target/plaac_headers.txt > target/_plaac_headers.haml

echo "ok (see target/_plaac_headers.haml)"

if ! type "dot" > /dev/null; then
  echo "dot not found, skipping genereration of hmm_default.png" 
else
  dot target/hmm_default.txt -Tpng > target/hmm_default.png
  echo "ok (see target/hmm_default.png)" 
fi


SCRIPT_DIR=`dirname $0`
cd $SCRIPT_DIR

echo cleanup
rm -rf target
mkdir target

echo build
javac -d target src/*.java

echo package
cp src/mainClass target/.
cp -r src/data target/.
cd target
jar cmf mainClass spewey.jar *
chmod +x spewey.jar

echo "ok (see target/spewey.jar)"

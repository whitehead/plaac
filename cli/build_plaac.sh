
SCRIPT_DIR=`dirname $0`
cd $SCRIPT_DIR

echo cleanup
rm -rf target
mkdir target

echo build
javac -d target src/*.java

echo package
cp src/mainClass target/.
cp -r src/util target/.
cd target
jar cmf mainClass plaac.jar *
chmod +x plaac.jar

echo "ok (see target/plaac.jar)"

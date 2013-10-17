
SCRIPT_DIR=`dirname $0`
cd $SCRIPT_DIR

echo cleanup
rm -rf target
mkdir target

echo build
javac -d target plaac/*.java

echo package
cp plaac/mainClass target/.
cp -r plaac/util target/.
cd target
jar cmf mainClass plaac.jar *
chmod +x plaac.jar

echo "ok (see target/plaac.jar)"

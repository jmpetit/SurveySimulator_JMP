#! /bin/bash

# Determine directory name
ID=`date +"%Y%m%d"`
d="../SurveySim${ID}"
da=`date +"%F"`
dt=`LC_ALL=us date +"%d %b %Y"`
df=`LC_ALL=us date +"%b %d/%Y"`
v="2.0"
intended_audience="OSSOS team"
#intended_audience="General Public"

# Create a clean distribution directory
if [ -d $d ]; then
    \rm -rf $d
fi
mkdir $d
if [ -e ../CurrentDistrib ]; then
    \rm ../CurrentDistrib
fi
if [ -h ../CurrentDistrib ]; then
    \rm ../CurrentDistrib
fi
ln -s ${d/./} ../CurrentDistrib
curdir=$(pwd)

# Copy fix files
cat > $d/README.version <<EOF

## Survey Simulator for OSSOSv12

## Survey simulator as of $dt

EOF
head -2 README.md > $d/README.md
cat >> $d/README.md <<EOF
$df release to $intended_audience
EOF
tail --line=+3 README.md >> $d/README.md
cp -a eupl* fortran python $d/
for s in CFEPS OSSOS All_r_Surveys All_Surveys; do
    mkdir -p $d/Surveys/$s
    cp Surveys/$s/* $d/Surveys/$s/
done
for s in CFEPS OSSOS All_r_Surveys; do
    \rm -f $d/Surveys/$s/README.formats
    ln -s ../All_Surveys/README.formats $d/Surveys/$s/README.formats
done
\rm -f $d/Simulator/F95/fortran/test*f95

cd ${curdir}

# Create the tarball
tar czf ../SurveySimulator-${ID}.tgz $d

exit

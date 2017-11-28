#$ -S /bin/bash
#$ -q cgc.q
#$ -cwd
#$ -V

umask 0007
umask
set

TIMESTAMP=`date '+%F.%H-%M-%S'`

echo this is a test at $TIMESTAMP >> /isilon/ddl/SS_0500181/foo.test.txt

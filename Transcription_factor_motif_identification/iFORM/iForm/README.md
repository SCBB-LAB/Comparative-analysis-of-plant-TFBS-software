Requirements：The GSL library is required，you need to install GSL library first。You can download from the follow url： http://ftp.kaist.ac.kr/gnu/gnu/gsl/ 
The latest version should work and version 2.2.1 has been tested. It is mentioned that the gsl library is install to /usr/local/lib by default, you should specify your valid path and add it to the LD_LIBRARY_PATH first.


INSTALL : after install the gsl library, use commandline $make to make the file, use $./iForm to excute the program
USAGE: iForm [options] <motif file> <sequence file>

More information is avaliable by using parameter -h
Troubleshooting:
In some situation the command "chmod -R 777 iForm" may need to be excuted after the make command.

#!/bin/csh
#
#  Purpose:
#
#    Create a GZIP'ed TAR file of the m_src/hammersley files.
#
#  Modified:
#
#    02 January 2006
#
#  Author:
#
#    John Burkardt
#
#  Move to the directory just above the "hammersley" directory.
#
cd $HOME/public_html/m_src
#
#  Delete any TAR or GZ file in the directory.
#
echo "Remove TAR and GZ files."
rm hammersley/*.tar
rm hammersley/*.gz
#
#  Create a TAR file of the "hammersley" directory.
#
echo "Create TAR file."
tar cvf hammersley_m_src.tar hammersley/*
#
#  Compress the file.
#
echo "Compress the TAR file."
gzip hammersley_m_src.tar
#
#  Move the compressed file into the "hammersley" directory.
#
echo "Move the compressed file into the directory."
mv hammersley_m_src.tar.gz hammersley
#
#  Say goodnight.
#
echo "The hammersley_m_src gzip file has been created."

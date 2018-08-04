#!/bin/sh
eval "$(ssh-agent -s)"
ssh-add
tmpdir=`mktemp -d`
if [ $# -gt 0 ]; then
  destdir=$1
else
  destdir=$HOME
fi
trap ":" INT

seconds=10 # or however many seconds you like

while sleep $seconds ; do
  rsync -c --progress --temp-dir="$tmpdir/" observer@epsilon:/home/observer/.jastrocam3/preview.fits $destdir
done
\rm -fr $tmpdir

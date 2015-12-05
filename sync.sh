#!/bin/sh

while true ; do
sshpass -p VW86Jetta rsync -c --temp-dir="/home/pi/Temp/rsync/data/temp/" zielona@beta:/home/zielona/.jastrocam3/preview.fits ~/Temp/rsync/data/
sshpass -p VW86Jetta rsync -c --temp-dir="/home/pi/Temp/rsync/data/temp/" zielona@beta:/home/zielona/dane/20150412/*.fits ~/Temp/rsync/data/
sleep 1 # or however many seconds you like
done

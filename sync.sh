#!/bin/sh

while true ; do
sshpass -p <pass> rsync -c --temp-dir="/home/pi/Temp/rsync/data/temp/" zielona@beta:/home/zielona/.jastrocam3/preview.fits ~/Temp/rsync/data/
#sshpass -p <PASSWORD> -c --temp-dir="/home/pi/Temp/rsync/data/temp/" zielona@beta:/home/zielona/dane/20151206/*.fits ~/Temp/rsync/data/
sleep 1 
done

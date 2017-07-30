#!/bin/sh
eval "$(ssh-agent -s)"
ssh-add
while true ; do
rsync -c --temp-dir="/home/ogloza/Programs/ObsAssistant/temp/temp/" observer@epsilon:/home/observer/.jastrocam3/preview.fits /home/ogloza/Programs/ObsAssistant/temp/data/
sleep 1 # or however many seconds you like
done

import sys
import glob
import shutil
import os
tempDir = '/home/pi/Temp/rsync/data/temp/'
dataDir = '/home/pi/Temp/rsync/data/'
ext = '*.temp'
os.system("sshpass -p VW86Jetta rsync zielona@beta:/home/zielona/.jastrocam3/preview.fits ~/Temp/rsync/data/temp/")
os.system("sshpass -p VW86Jetta rsync zielona@beta:/home/zielona/dane/20150412/*.fits ~/Temp/rsync/data/temp/")
tempList = sorted(glob.glob(tempDir + ext))
dataList = sorted(glob.glob(dataDir + ext))

while True:
    diff = set(tempList).intersection(dataList)
    if len(diff) > 0:
        for fi in diff:
            shutil.copy(tempDir + fi, dataDir)

    os.system("sshpass -p VW86Jetta rsync zielona@beta:/home/zielona/.jastrocam3/preview.fits ~/Temp/rsync/data/temp/")
    os.system("sshpass -p VW86Jetta rsync zielona@beta:/home/zielona/dane/20150412/*.fits ~/Temp/rsync/data/temp/")
    tempList = sorted(glob.glob(tempDir + ext))
    dataList = sorted(glob.glob(dataDir + ext))



# to start the container
docker run -it --platform linux/amd64 --rm \
    -v /Users/josephchen/Library/CloudStorage/Box-Box/josephchen/MnM/MnMscripts_forpub/data/MnM_dFC:/out \
    -v /Users/josephchen/Library/CloudStorage/Box-Box/josephchen/MnM/MnMscripts_forpub/data/events:/ratings \
    -v /Users/josephchen/data/MnM/MnM-XCPD-NoGSRNoCensor:/data \
    -v /Users/josephchen/Library/CloudStorage/Box-Box/josephchen/MnM/MnMscripts_forpub/methods1_TVDN:/script \
    huaqingjin/tvdn /bin/bash

# navigate to output folder
cd /out/TVDN
git clone https://github.com/JINhuaqing/TVDN.git
cd TVDN
cp /script/dFC_TVDNDocker_2BroadAnalysis.py .
cp /script/dFC_TVDNDocker_3Summarise_masterdf.py .
cp /script/dFC_TVDNDocker_4ActualAnalysis_Graph.py .

# or run the code
python dFC_TVDNDocker_2BroadAnalysis.py
python dFC_TVDNCovker_3Summarise_masterdf.py
python dFC_TVDNDocker_4ActualAnalysis_Graph.py


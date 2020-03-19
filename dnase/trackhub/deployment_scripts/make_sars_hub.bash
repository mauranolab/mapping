#!/bin/bash

#####################################################################
# Clean out old hub:

rm -rf /home/cadlej01/public_html/trackhub_sars
mkdir /home/cadlej01/public_html/trackhub_sars

#####################################################################

src_dev/update_tracks.bash SARS /home/cadlej01/public_html/trackhub_sars "SARS" "SARS Hub" \
                      https://mauranolab:chromatin@cascade.isg.med.nyu.edu/~cadlej01

# Hub address will be:  https://mauranolab:chromatin@cascade.isg.med.nyu.edu/~cadlej01/trackhub_sars/hub.txt


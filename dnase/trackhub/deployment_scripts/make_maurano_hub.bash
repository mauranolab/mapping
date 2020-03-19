#!/bin/bash

rm -rf /home/cadlej01/public_html/trackhub
mkdir /home/cadlej01/public_html/trackhub

# src_prod/update_tracks.bash MAURANOLAB /home/cadlej01/public_html/trackhub "Maurano Lab" "Maurano Lab Hub" \
#                       https://mauranolab:chromatin@cascade.isg.med.nyu.edu/~cadlej01

src_dev/update_tracks.bash MAURANOLAB /home/cadlej01/public_html/trackhub "Maurano Lab" "Maurano Lab Hub" \
                      https://mauranolab:chromatin@cascade.isg.med.nyu.edu/~cadlej01

# Hub address will be:  https://mauranolab:chromatin@cascade.isg.med.nyu.edu/~cadlej01/trackhub/hub.txt


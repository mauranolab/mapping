#!/bin/bash

#####################################################################
# Clean out old hub:

rm -rf /vol/cegs/public_html/trackhub/hg38
rm -rf /vol/cegs/public_html/trackhub/mm10
rm -rf /vol/cegs/public_html/trackhub/rn6
rm -rf /vol/cegs/public_html/trackhub/cegsvectors

rm -f  /vol/cegs/public_html/trackhub/genomes.txt
rm -f  /vol/cegs/public_html/trackhub/description.html
rm -f  /vol/cegs/public_html/trackhub/hub.txt
rm -f  /vol/cegs/public_html/trackhub/README

rm -f  /vol/cegs/public_html/trackhub/publicdata
rm -f  /vol/cegs/public_html/trackhub/aggregations
rm -f  /vol/cegs/public_html/trackhub/mapped

#####################################################################

src_prod/update_tracks.bash CEGS /vol/cegs/public_html/trackhub "CEGS" "CEGS Hub" \
                                 https://cegs:bigdna@cascade.isg.med.nyu.edu/cegs

# Hub address will be:  https://cegs:bigdna@cascade.isg.med.nyu.edu/cegs/trackhub/hub.txt


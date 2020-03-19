#!/bin/bash

#####################################################################
# Clean out old hub:

rm -rf /vol/cegs/public_html/trackhub_dev/hg38
rm -rf /vol/cegs/public_html/trackhub_dev/mm10
rm -rf /vol/cegs/public_html/trackhub_dev/rn6
rm -rf /vol/cegs/public_html/trackhub_dev/cegsvectors

rm -f  /vol/cegs/public_html/trackhub_dev/genomes.txt
rm -f  /vol/cegs/public_html/trackhub_dev/description.html
rm -f  /vol/cegs/public_html/trackhub_dev/hub.txt
rm -f  /vol/cegs/public_html/trackhub_dev/README

rm -f  /vol/cegs/public_html/trackhub_dev/publicdata
rm -f  /vol/cegs/public_html/trackhub_dev/aggregations
rm -f  /vol/cegs/public_html/trackhub_dev/mapped

#####################################################################

src_dev/update_tracks.bash CEGS /vol/cegs/public_html/trackhub_dev "CEGS Dev" "CEGS Development Hub" \
                                https://cegs:bigdna@cascade.isg.med.nyu.edu/cegs

# Hub address will be:  https://cegs:bigdna@cascade.isg.med.nyu.edu/cegs/trackhub_dev/hub.txt


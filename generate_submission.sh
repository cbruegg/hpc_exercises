#!/bin/bash

if [[ "$#" -ne 1 ]]; then
    echo "You need to specify the lab number!"
    exit 1
fi

NUM=$1
PAD_NUM=`printf "%02d" ${NUM}`

ADDITIONAL_FILES=""
if (($NUM == 6)); then
    ADDITIONAL_FILES="${ADDITIONAL_FILES} input*.ppm"
fi

TRACKED=`git ls-tree -r HEAD --name-only | awk 1 ORS=' ' | tr " " "\n" | sort -g`
LAB_FILES=`echo "ManualMakefile lab${PAD_NUM}* lab${PAD_NUM}*/* ${ADDITIONAL_FILES}" | tr " " "\n" | sort -g`
INTERSECTION_FILES=`comm -12 <(echo ${TRACKED} | tr " " "\n") <(echo ${LAB_FILES} | tr " " "\n")`

zip submission_${PAD_NUM}.zip ${INTERSECTION_FILES}
#!/bin/bash

NAME=allelecount_ggplot
VERSION=0.1

docker build -t ${NAME}:${VERSION} .

singularity build ${NAME}.${VERSION}.sif docker-daemon://${NAME}:${VERSION}

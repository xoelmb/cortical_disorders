#!/bin/sh

conda create -y -p ./conda_env "r-base>=4.1" r-essentials r-arrow cmake libgit2 r-systemfonts r-textshaping r-ragg

conda activate ./conda_env

Rscript -e 'source("init.R", echo=TRUE)'
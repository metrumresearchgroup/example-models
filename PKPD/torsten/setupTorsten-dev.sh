# setupTorsten-dev
# v 1.3
# Run this bash file in the directory where you want to install cmdstan
# with the Torsten functions. 
#
# update 1.3: download cmdStan v2.16.0.
# update dev: download dev version of Torsten instead of last release.
# update 1.2: download cmdStan v2.14.0 and the master version of 
#             torsten-stan and torsten-math.
# update 1.1: download last version of cmdStan/dev with which Torsten
#             was tested.

#!/bin/bash
git clone https://github.com/stan-dev/cmdstan.git
cd cmdstan
git checkout release/v2.16.0
git clone https://github.com/charlesm93/stan.git
# rm -r stan_2.16.0
# mv stan stan_2.16.0
# cd stan_2.16.0
cd stan
git checkout torsten-develop
# mkdir lib
cd lib
git rm -r stan_math
git clone https://github.com/charlesm93/math.git
mv math stan_math
cd stan_math
git checkout torsten-develop

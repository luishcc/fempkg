#!/bin/bash

dir=$HOME/fempkg


sim_dir=$dir/sim
fem_env=$dir/fem.yml



if ! test -d "$sim_dir/results"; then
  mkdir $sim_dir/results
fi

ans="n"
if ! test -d "$HOME/miniconda3"; then
  printf "\nInstall Miniconda? y/n\n"
  read ans
else
  printf "\nMiniconda already installed\n"
fi

if [[ "$ans" = "y" ]]; then
  printf  "\nThe miniconda Installation should be made in the $HOME directory
  and say \"yes\"  to the option of running conda init \n
  Download the installer? y/n"
  read ans2
  if [[ "$ans2" = "y" ]]; then
    wget -P $HOME https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
  fi
  chmod +x $HOME/Miniconda3-latest-Linux-x86_64.sh
  $HOME/Miniconda3-latest-Linux-x86_64.sh
  conda config --set auto_activate_base false
  printf "\nInstallation Done\nSetting fem.yml\n"
fi

source $HOME/miniconda3/etc/profile.d/conda.sh
conda env create --file $fem_env

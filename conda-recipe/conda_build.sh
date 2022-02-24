PKG_NAME=mtsv-tools
USER=tfursten

OS=linux-64
mkdir ~/conda-bld
conda config --set anaconda_upload no
export CONDA_BLD_PATH=~/conda-bld
export VERSION="2.0.0"
conda build ..
conda install mtsv_tools --use-local

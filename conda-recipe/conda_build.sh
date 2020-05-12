PKG_NAME=mtsv-tools
USER=tfursten

OS=linux-64
mkdir ~/conda-bld
conda config --set anaconda_upload no
export CONDA_BLD_PATH=~/conda-bld
export VERSION="1.0.1"
conda build ..
conda install your-package --use-local
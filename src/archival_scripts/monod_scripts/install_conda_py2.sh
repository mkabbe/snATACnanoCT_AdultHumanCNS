set -e  # Stop on error

CONDA_ENV_PY2=encode-atac-seq-pipeline-python2
REQ_TXT_PY2=requirements_py2.txt

conda --version  # check if conda exists

echo "=== Installing pipeline's Conda environments ==="
conda create -n ${CONDA_ENV_PY2} --file ${REQ_TXT_PY2} -y -c defaults -c r -c bioconda -c conda-forge


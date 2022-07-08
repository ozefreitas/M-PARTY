FROM continuumio/miniconda3

RUN conda update -n base -c defaults -y conda \
&& git clone https://github.com/ozefreitas/PlastEDMA.git \
&& conda install -c conda-forge -y mamba \
&& mamba env update -f PlastEDMA/ci/ci_environment.yml --name base \
&& bash PlastEDMA/ci/ci_build.sh \
&& conda clean --all -y

CMD [ "python", "bin/workflow/plastedma.py" ]

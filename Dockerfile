FROM continuumio/miniconda3

RUN git clone https://github.com/pg42872/PlastEDMA.git \
&& conda install -c conda-forge -y mamba \
&& mamba env update --file PlastEDMA/workflow/envs/environment.yml --name base \
&& bash PlastEDMA/ci/ci_build.sh

CMD [ "python", "bin/plastedma.py" ]

FROM bentoml/model-server:0.11.0-py311
MAINTAINER ersilia

RUN conda install -c conda-forge libstdcxx-ng
RUN pip install dilipred==4.0.9

WORKDIR /repo
COPY . /repo

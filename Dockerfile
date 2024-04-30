FROM bentoml/model-server:0.11.0-py310
MAINTAINER ersilia

RUN pip install rdkit-pypi==2022.9.5
RUN pip install protobuf==3.18.3
RUN pip install tensorboardx==2.0
RUN pip install chemprop==1.5.2
RUN pip install tensorboard==2.11.0

# Clean up the temporary directory
RUN rm -rf /temp_repo

# Clone the repository and checkout the desired commit
RUN mkdir /temp_repo && cd /temp_repo && git clone https://github.com/MolecularAI/REINVENT4 . && rm -rf priors && git add . && git commit -m "remove prior" && git checkout ee87830e800af22f17b05027bb9378950d5980f7

# Install the package using pip
RUN pip install /temp_repo --extra-index-url=https://pypi.anaconda.org/OpenEye/simple --extra-index-url https://download.pytorch.org/whl/cu113

WORKDIR /repo
COPY . /repo

FROM bentoml/model-server:0.11.0-py310
MAINTAINER ersilia

RUN python -m pip install --upgrade pip

RUN pip install torch==2.9.1 --index-url https://download.pytorch.org/whl/cpu
RUN pip install torchvision==0.24.1 --index-url https://download.pytorch.org/whl/cpu

RUN pip install chemprop==1.5.2 --extra-index-url https://pypi.anaconda.org/OpenEye/simple
RUN pip install descriptastorus==2.8.0
RUN pip install funcy==2.0
RUN pip install matplotlib==3.8.4
RUN pip install mmpdb==2.1
RUN pip install MolVS==0.1.1
RUN pip install numpy==1.26.4
RUN pip install pandas==2.2.2
RUN pip install polars==1.37.1
RUN pip install Pillow==10.4.0
RUN pip install pumas==1.3.0
RUN pip install pydantic==2.11.0
RUN pip install pytest==8.2.2
RUN pip install pytest-mock==3.14.0
RUN pip install python-dotenv==1.0.1
RUN pip install PyYAML==6.0.2
RUN pip install rdkit==2025.9.2
RUN pip install requests==2.32.3
RUN pip install requests_mock==1.12.1
RUN pip install tenacity==8.5.0
RUN pip install tensorboard==2.20.0
RUN pip install tomli==2.4.0
RUN pip install tqdm==4.66.5
RUN pip install apted==1.0.3
RUN pip install typing_extensions==4.13.2
RUN pip install xxhash==3.5.0
RUN pip install pathos==0.3.2
RUN pip install git+https://github.com/MolecularAI/REINVENT4@v4.2.6 --no-deps --extra-index-url https://pypi.anaconda.org/OpenEye/simple

WORKDIR /repo
COPY . /repo

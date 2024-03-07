# Use Ubuntu 22.04 as the base image
FROM ubuntu:22.04

# Install sudo
RUN apt-get update && \
    apt-get install -y sudo && \
    rm -rf /var/lib/apt/lists/*

# Install system dependencies
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    python3.10 \
    python3-venv \
    python3-pip \
    python3-setuptools \
    python3-wheel \
    python3-dev \
    libqt5multimedia5-plugins \
    build-essential \
    gcc \
    pyqt5-dev \
    xvfb \
    && rm -rf /var/lib/apt/lists/*

# Set the default Python version to use
RUN ln -s /usr/bin/python3.10 /usr/bin/python

# # Set environment variables
ENV QT_DEBUG_PLUGINS=1
ENV DISPLAY=:1

# Create and set working directory
WORKDIR /go-rxr

## Copy project files into the container
# Folders
COPY ./DATA /go-rxr/DATA
COPY ./DOCS /go-rxr/DOCS
COPY ./SCRIPTS /go-rxr/SCRIPTS
COPY ./TESTS /go-rxr/TESTS
COPY ./UTILS /go-rxr/UTILS
# Files
COPY ./requirements.txt /go-rxr/requirements.txt
COPY ./setup.py /go-rxr/setup.py
COPY ./setup_reflectivity.py /go-rxr/setup_reflectivity.py
COPY ./GUI_GO.py /go-rxr/GUI_GO.py
COPY ./Pythonreflectivity.cpp /go-rxr/Pythonreflectivity.cpp
COPY ./Pythonreflectivity.pyx /go-rxr/Pythonreflectivity.pyx
COPY ./form_factor.pkl /go-rxr/form_factor.pkl
COPY ./form_factor_magnetic.pkl /go-rxr/form_factor_magnetic.pkl
COPY ./Ti34OpsPython.pkl /go-rxr/Ti34OpsPython.pkl

# # Create and activate virtual environment
# RUN python3.10 -m venv venv-go-rxr
# RUN /bin/bash -c "source venv-go-rxr/bin/activate"

# Install required Python libraries
RUN python3.10 -m venv venv-go-rxr && \
    . venv-go-rxr/bin/activate && \
    pip install -r requirements.txt

# Install required Python libraries
RUN pip install -r requirements.txt

# Install Python reflectivity
RUN python setup_reflectivity.py install

# # Set default command to run the application
# CMD ["python", "GUI_GO.py"]
# CMD [ "ls -l" ]
# Use Ubuntu 22.04 as the base image
FROM ubuntu:22.04

# Install system dependencies
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    python3.10 \
    python3-venv \
    python3-pip \
    python3-setuptools \
    python3-wheel \
    libqt5multimedia5-plugins \
    && rm -rf /var/lib/apt/lists/*

# Set the default Python version to use
RUN ln -s /usr/bin/python3.10 /usr/bin/python

# Create and set working directory
WORKDIR /go-rxr

# Copy project files into the container
# COPY ../requirements.txt /go-rxr/requirements.txt
# COPY ../

# COPY ../DATA /go-rxr/DATA
# COPY ../DOCS /go-rxr/DOCS
# COPY ../SCRIPTS /go-rxr/SCRIPTS
# COPY ../TESTS /go-rxr/TESTS
# COPY ../UTILS /go-rxr/UTILS
COPY ../. /go-rxr

# # Create and activate virtual environment
# RUN python3.10 -m venv venv-go-rxr
# RUN /bin/bash -c "source venv-go-rxr/bin/activate"

# # Install Python libraries
# RUN python setup.py install

# # Install additional Python libraries
# RUN pip install Pillow

# # Install Python reflectivity
# RUN python setup_reflectivity.py install

# # Expose any necessary ports
# # EXPOSE <port_number>

# # Set default command to run the application
# CMD ["python", "GUI_GO.py"]
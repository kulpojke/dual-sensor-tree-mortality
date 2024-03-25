FROM continuumio/miniconda3

# Set the working directory in the container
WORKDIR /app

# Copy the environment.yml file into the container at /app
COPY environment.yml /app/environment.yml

# Set up the Conda environment
RUN echo "source activate $(head -1 /app/environment.yml | cut -d' ' -f2)" > ~/.bashrc
ENV PATH /opt/conda/bin:$PATH

# Create the Conda environment using the environment.yml file
RUN conda env create -f /app/environment.yml

# Activate the environment
SHELL ["conda", "run", "-n", "$(head -1 /app/environment.yml | cut -d' ' -f2)", "/bin/bash", "-c"]

# The code or application to run in the container should be added here
# COPY your-application /app

# The default command to run when starting the container
# CMD ["your-application-command"]

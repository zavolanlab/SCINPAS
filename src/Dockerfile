


# build our own image based on base image 
# FROM gcc:11.1.0 AS dependencies
# FROM python:3.9.12 AS base
FROM ubuntu:latest

# set the working directory in the container. This is about inside the container (inside your local repository).
WORKDIR /app

# Copy copy the dependencies file to destination within docker image or container itself. (In this case, our working directory)
COPY requirements.txt /app

# install gcc-11
# RUN add-apt-repository ppa:ubuntu-toolchain-r/test
RUN apt-get update \
    && apt-get install -y gcc-11 g++-11

# install python and pip
# but why not specific version like python3.9.12?
RUN apt-get update && apt-get install -y \
    python3 \
    python3-pip

# install dependencies
RUN pip3 install -r requirements.txt


# move python script
# At this point, we have an image that is based on Python version 3.9.12 and we have installed our dependencies. The next step is to add our source code into the image.
# This COPY command (COPY . .) takes all the files located in the current directory and copies them into the image. 
# This COPY command (COPY src/ ) copies all content of the local src directory to the working directory
# COPY src/ .

# Now, all we have to do is to tell Docker what command we want to run when our image is executed inside a container. We do this using the CMD command. 
# command to run on container start

# CMD ["python", "./get_barchart_full.py"]
############################################################################################################################################################################

# to build a container (from image -> container):
# 1) go to folder where docker file is
# 2) docker build . (. is current path)

# to see if image is build use this command: dokcer images

# if you want to build container with a tag use this command: docker build -tag name_of_the_image
# use this command docker build -t zavolab/polcat_python:v1.0.0 .
# the last dot indicates that you are doing this in current directory. (you can give absolute path here instead of dot)

# A Docker container is a virtualized runtime environment that provides isolation capabilities for separating the execution of applications from the underpinning system. 
# It’s an instance of a Docker image.

# Docker defines several container states, such as created, restarting, running, paused, exited, and dead. 
# Since several states are possible, and a container is just an instance of the image, a container does not need to be running.

# Every time Docker creates a container from an image, it places a thin read-write layer on top of the image. 
# This writable layer allows for changes to be made to the container, 
# as the lower layers in the image are unchangeable. It also stores any changes made to the container during its entire runtime.


FROM python:3.13.3-bullseye

# Running apt update to ensure the package list is up to date
RUN apt-get update && apt-get install -y --no-install-recommends graphviz

# Create directory for app
RUN mkdir -p /usr/src/app
WORKDIR /usr/src/app

# Configure Python environment
COPY requirements.txt /usr/src/app/
RUN pip install --no-cache-dir -r requirements.txt
COPY . /usr/src/app
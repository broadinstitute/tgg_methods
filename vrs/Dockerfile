# Dockerfile for creating an image that can run the VRS Annotation Script inside of it
# includes libpqdev and pythondev, plus psycopg2, ga4gh.vrs[extras] , and seqrepo

FROM 'google/cloud-sdk:slim'

RUN apt-get -y update && apt-get install -y --no-install-recommends \
    apt-utils \
    docker-compose \
    g++ \
    gawk \
    less \
    libbz2-dev \
    libcurl4-openssl-dev \
    liblzma-dev \
    libncurses5-dev \
    man-db \
    pkg-config \
    python3-venv \
    software-properties-common \
    unzip \
    wget \
    zlib1g-dev \
    libpq-dev \ 
    python-dev\
    && \
	# clean up apt cache
    rm -rf /var/lib/apt/lists/*

# Install Java 8 for hail
RUN wget -qO - https://adoptopenjdk.jfrog.io/adoptopenjdk/api/gpg/key/public | apt-key add - \ 
    && add-apt-repository -y https://adoptopenjdk.jfrog.io/adoptopenjdk/deb/ \ 
    && apt-get update && apt-get install -y adoptopenjdk-8-hotspot

RUN mkdir tools 
WORKDIR /tools



# Install vcftools
ENV VCFTOOLS_VERSION="0.1.16"
RUN wget https://github.com/vcftools/vcftools/releases/download/v${VCFTOOLS_VERSION}/vcftools-${VCFTOOLS_VERSION}.tar.gz \
    && tar xf vcftools-${VCFTOOLS_VERSION}.tar.gz \
    && rm vcftools-${VCFTOOLS_VERSION}.tar.gz \
	&& cd vcftools-${VCFTOOLS_VERSION} \
	&& ./configure \
	&& make \
	&& make install \
	&& make clean


# Install python packages
RUN apt-get update \
	&& apt-get dist-upgrade -y \
	&& apt-get install -y --no-install-recommends\
    libc6-dev \
    libffi-dev \
    libgdbm-dev \
    liblapack-dev \
    liblapack3 \ 
    libncursesw5-dev \
    libopenblas-base \
    libopenblas-dev \
    libsqlite3-dev \
    libssl-dev \
    openssl \
    python3-smbus \
    python3 python3-pip

# Upgrade pip to latest version
RUN python3 -m pip install --upgrade pip

# Install hail and other python libraries
# Including packages needed for VRS annotation functionality
# Also include version for ga4gh.vrs-python to be most recent and supporting their new variant schema  
# Please update the constant VRS_VERSION in vrs_annotation_batch.py accordingly when GA4GH_VRS_VERSION is changed
ENV HAIL_VERSION="0.2.113" 
ENV GA4GH_VRS_VERSION="0.8.4"
RUN python3 --version
RUN python3 -m pip install \
    wheel \
    pypandoc \
    hail==${HAIL_VERSION} \
    scipy \
    numpy \
    pandas \
    matplotlib \
    seaborn \
    ipython \
    pybedtools \
    dill \
    psycopg2 \
    ga4gh.vrs[extras]==${GA4GH_VRS_VERSION}  \
    seqrepo \
    click

# Install GCS Connector
RUN export SPARK_HOME=$(find_spark_home.py) && \
    curl https://storage.googleapis.com/hadoop-lib/gcs/gcs-connector-hadoop2-2.0.1.jar \
         >$SPARK_HOME/jars/gcs-connector-hadoop2-2.0.1.jar && \
    mkdir -p $SPARK_HOME/conf && \
    touch $SPARK_HOME/conf/spark-defaults.conf && \
    sed -i $SPARK_HOME/conf/spark-defaults.conf \
        -e 's:spark\.hadoop\.google\.cloud\.auth\.service\.account\.enable.*:spark.hadoop.google.cloud.auth.service.account.enable true:' \
        -e 's:spark\.hadoop\.google\.cloud\.auth\.service\.account\.json\.keyfile.*:spark\.hadoop\.google\.cloud\.auth\.service\.account\.json\.keyfile /gsa-key/key.json:'

# Copy in SeqRepo from local file structure
COPY /usr/local/share/seqrepo/2018-11-26/ /tmp/local-seqrepo

WORKDIR /home

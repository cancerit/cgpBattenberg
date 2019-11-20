FROM  quay.io/wtsicgp/dockstore-cgpmap:3.1.4 as builder

USER  root

# ALL tool versions used by opt-build.sh
ENV VER_ALLELECOUNT="4.1.0"
ENV VER_CGPVCF="v2.2.1"
ENV VER_IMPUTE2="v2.3.0"
ENV VER_VCFTOOLS="0.1.16"

RUN apt-get -yq update
RUN apt-get install -yq --no-install-recommends \
locales \
g++ \
make \
gcc \
pkg-config \
python \
zlib1g-dev \
r-base \
libbz2-dev \
liblzma-dev \
libcurl4-openssl-dev \
libssl-dev \
libxml2-dev \
libssh2-1-dev

RUN locale-gen en_US.UTF-8
RUN update-locale LANG=en_US.UTF-8

ENV OPT /opt/wtsi-cgp
ENV PATH $OPT/bin:$OPT/biobambam2/bin:$PATH
ENV PERL5LIB $OPT/lib/perl5
ENV LD_LIBRARY_PATH $OPT/lib
ENV LC_ALL en_US.UTF-8
ENV LANG en_US.UTF-8
ENV R_LIBS $OPT/R-lib
ENV R_LIBS_USER $R_LIBS
ENV R_PROFILE_USER $OPT/config/Rprofile

COPY build/Rprofile $OPT/config/Rprofile

# build tools from other repos
ADD build/opt-build.sh build/
RUN bash build/opt-build.sh $OPT

# build the tools in this repo, separate to reduce build time on errors
COPY . .
RUN bash build/opt-build-local.sh $OPT

FROM ubuntu:16.04

LABEL maintainer="cgphelp@sanger.ac.uk" \
      uk.ac.sanger.cgp="Cancer, Ageing and Somatic Mutation, Wellcome Trust Sanger Institute" \
      version="3.5.1" \
      description="cgpCaVEManWrapper docker"

RUN apt-get -yq update
RUN apt-get install -yq --no-install-recommends \
apt-transport-https \
locales \
curl \
ca-certificates \
libperlio-gzip-perl \
bzip2 \
psmisc \
time \
zlib1g \
liblzma5 \
libncurses5 \
p11-kit \
r-base \
unattended-upgrades && \
unattended-upgrade -d -v && \
apt-get remove -yq unattended-upgrades && \
apt-get autoremove -yq

RUN locale-gen en_US.UTF-8
RUN update-locale LANG=en_US.UTF-8

ENV OPT /opt/wtsi-cgp
ENV PATH $OPT/bin:$OPT/biobambam2/bin:$PATH
ENV PERL5LIB $OPT/lib/perl5
ENV LD_LIBRARY_PATH $OPT/lib
ENV LC_ALL en_US.UTF-8
ENV LANG en_US.UTF-8
ENV R_LIBS $OPT/R-lib
ENV R_LIBS_USER $R_LIBS
ENV R_PROFILE_USER $OPT/config/Rprofile

RUN mkdir -p $OPT
COPY --from=builder $OPT $OPT

## USER CONFIGURATION
RUN adduser --disabled-password --gecos '' ubuntu && chsh -s /bin/bash && mkdir -p /home/ubuntu

USER    ubuntu
WORKDIR /home/ubuntu

CMD ["/bin/bash"]

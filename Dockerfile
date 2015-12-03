FROM rocker/r-base:latest

WORKDIR /grip-r

RUN apt-get update && apt-get install -y openssl libcurl4-openssl-dev curl libxml2-dev libssl-dev

ADD ./inst/bash/install-package-dependencies.sh /grip-r/inst/bash/install-package-dependencies.sh

RUN ./inst/bash/install-package-dependencies.sh

ADD ./ /grip-r

RUN R --no-save --quiet -e 'devtools::document()'
RUN R CMD INSTALL --no-multiarch --with-keep.source /grip-r
RUN R CMD build /grip-r

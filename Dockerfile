FROM rocker/r-ver:4.0.3

RUN apt-get update && \
  apt-get install -y libxml2-dev libcurl4-openssl-dev libssl-dev
RUN R -e "install.packages('devtools')"
RUN R -e "install.packages('collections')"
RUN R -e "devtools::install_github('ctlab/linseed')"
RUN R -e "install.packages('plotly')"
RUN R -e "install.packages('matlib')"
RUN R -e "install.packages('matrixcalc')"
RUN R -e "install.packages('optparse')"
RUN R -e "install.packages('yaml')"
RUN R -e "install.packages('nnls')"
RUN R -e "install.packages('Rcpp')"
RUN R -e "install.packages('RcppArmadillo')"


WORKDIR /app
COPY ./SinkhornNNLSLinseedC.R ./SinkhornNNLSLinseedC.R
COPY ./LinseedMetadata.R ./LinseedMetadata.R
COPY ./RunDeconvolution.R ./RunDeconvolution.R
COPY ./pipeline.cpp ./pipeline.cpp
CMD ["Rscript", "RunDeconvolution.R"]

# Container with R, Bioconductor and every package the code bank needs.
#
# Build:   docker build -t rnaseq-codebank .
# Test:    docker run --rm rnaseq-codebank Rscript tests/run_tests.R
# Run:     docker run --rm -v "$PWD/runs:/codebank/runs" -v "/path/to/data:/data" rnaseq-codebank \
#            Rscript scripts/run_pipeline.R --counts_file /data/counts.csv --metadata_file /data/samples.csv \
#            --condition condition --reference_level control --test_level treated
#
# The Bioconductor image already includes R, pandoc and system libraries.
# For exact reproducibility, replace "latest" with a fixed release tag
# (e.g. RELEASE_3_22; see https://hub.docker.com/r/bioconductor/bioconductor_docker/tags)
# and record it in CHANGELOG.md.
ARG BIOC_TAG=latest
FROM bioconductor/bioconductor_docker:${BIOC_TAG}

WORKDIR /codebank
COPY scripts/00_setup_environment.R scripts/00_setup_environment.R
RUN Rscript scripts/00_setup_environment.R

COPY . .
RUN Rscript scripts/make_example_data.R

CMD ["Rscript", "scripts/run_pipeline.R", "--help"]

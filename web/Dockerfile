# temporary container for building CLI plaac.jar

FROM docker.io/library/ruby:3.3-bookworm AS plaac-builder

ENV DEBIAN_FRONTEND=noninteractive

ARG PLAAC_VERSION

ENV SETUPTOOLS_SCM_PRETEND_VERSION_FOR_PLAAC=${PLAAC_VERSION}

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        openjdk-17-jdk-headless \
        python3 \
        python3-venv \
	graphviz \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /build

COPY cli/ ./cli/

RUN cd cli && ./build_plaac.sh

# main container for the web application

FROM docker.io/library/ruby:3.3-bookworm

LABEL org.opencontainers.image.source="https://github.com/whitehead/plaac"
LABEL org.opencontainers.image.description "PLAAC web server container image"
LABEL org.opencontainers.image.licenses="MIT"

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        openjdk-17-jre-headless \
	nodejs \
        r-base \
        build-essential \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app/web

COPY web/Gemfile web/Gemfile.lock ./

RUN bundle install

COPY web/ ./

# get the plaac.jar from the build container to include in the webserver

COPY --from=plaac-builder /build/cli/target/plaac.jar ./bin/plaac.jar
COPY --from=plaac-builder /build/cli/target/_plaac_headers.haml ./views/_plaac_headers.haml
COPY --from=plaac-builder /build/cli/target/hmm_default.png ./public/hmm_default.png

RUN mkdir -p tmp logs /tmp/clusty

EXPOSE 4567

HEALTHCHECK --interval=30s --timeout=5s --start-period=10s --retries=3 \
  CMD ruby -rnet/http -e 'exit(Net::HTTP.get_response(URI("http://127.0.0.1:4567/")).is_a?(Net::HTTPSuccess) ? 0 : 1)'

CMD ["bundle", "exec", "rackup", "config.ru", "-p", "4567", "-o", "0.0.0.0"]

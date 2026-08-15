FROM docker.io/library/ruby:3.3-bookworm

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        openjdk-17-jre-headless \
	nodejs \
        r-base \
        build-essential \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app

COPY Gemfile Gemfile.lock ./

RUN bundle install

COPY . .

RUN mkdir -p tmp logs /tmp/clusty

EXPOSE 4567

HEALTHCHECK --interval=30s --timeout=5s --start-period=10s --retries=3 \
  CMD ruby -rnet/http -e 'exit(Net::HTTP.get_response(URI("http://127.0.0.1:4567/")).is_a?(Net::HTTPSuccess) ? 0 : 1)'

CMD ["bundle", "exec", "rackup", "config.ru", "-p", "4567", "-o", "0.0.0.0"]
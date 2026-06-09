FROM debian:stable-slim

RUN apt-get update \
    && apt-get install -y --no-install-recommends build-essential zlib1g-dev ca-certificates \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /opt/tsoclip
COPY . .
RUN make clean && make && make check && make install PREFIX=/usr/local

ENTRYPOINT ["tsoclip"]
CMD ["--help"]

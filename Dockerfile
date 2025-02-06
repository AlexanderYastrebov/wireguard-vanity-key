FROM golang:latest AS builder
WORKDIR /app
COPY . .
RUN CGO_ENABLED=0 go build

FROM scratch
LABEL org.opencontainers.image.url=https://github.com/AlexanderYastrebov/wireguard-vanity-key \
  org.opencontainers.image.licenses=BSD-3-Clause
COPY --from=builder /app/wireguard-vanity-key /wireguard-vanity-key
ENTRYPOINT ["/wireguard-vanity-key"]
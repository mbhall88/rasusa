FROM alpine

ADD target/docker/rasusa /bin/

CMD ["/bin/rasusa"]

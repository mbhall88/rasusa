FROM scratch

ADD target/docker/rasusa /bin/

CMD ["/bin/rasusa"]
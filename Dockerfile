FROM alpine

RUN apk update

RUN apk upgrade

RUN apk add bash

ADD target/docker/rasusa /bin/

CMD ["/bin/rasusa"]

TAG1 = weisburd/image-name:v1
TAG2 = weisburd/image-name:latest

all: build push

build:
	docker build -t $(TAG1) -t $(TAG2) .

push:
	docker push $(TAG1)
	docker push $(TAG2)

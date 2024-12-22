### auto-run

### sc-recounter-download

Env vars

```bash
IMG_NAME=sc-recounter-run
IMG_VERSION=0.1.0
REGION="us-east1"
PROJECT="c-tc-429521"
```

Build

> from the base directory of the repository

```bash
docker build \
  --file docker/${IMG_NAME}/Dockerfile \
  --build-arg CONDA_ENV_YAML=docker/${IMG_NAME}/environment.yml \
  --platform linux/amd64 \
  --tag ${IMG_NAME}:${IMG_VERSION} \
  .
```

Run the image

```bash
docker run -it --rm \
  -u $(id -u):$(id -g) \
  -v ${PWD}:/data \
  -v ${HOME}/.gcp/:/.gcp \
  --env GOOGLE_APPLICATION_CREDENTIALS=/.gcp/c-tc-429521-6f6f5b8ccd93.json \
  --platform linux/amd64 \
  ${IMG_NAME}:${IMG_VERSION} \
  --help
```

Push

```bash
docker tag ${IMG_NAME}:${IMG_VERSION} \
  ${REGION}-docker.pkg.dev/${PROJECT}/${IMG_NAME}/${IMG_NAME}:${IMG_VERSION} \
  && docker push ${REGION}-docker.pkg.dev/${PROJECT}/${IMG_NAME}/${IMG_NAME}:${IMG_VERSION}
```



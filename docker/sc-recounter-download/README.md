### sc-recounter-download

Build

> from the base directory of the repository

```bash
IMG_NAME=sc-recounter-download
IMG_VERSION=0.1.0
docker build \
  --file docker/${IMG_NAME}/Dockerfile \
  --platform linux/amd64 \
  --tag ${IMG_NAME}:${IMG_VERSION} \
  .
```

Run

```bash
docker run -it --rm \
  -u $(id -u):$(id -g) \
  -v ${PWD}:/data \
  -v ${HOME}/.gcp/:/.gcp \
  --env GOOGLE_APPLICATION_CREDENTIALS="/.gcp/c-tc-429521-6f6f5b8ccd93.json" \
  --env GCP_PROJECT_ID="c-tc-429521" \
  --platform linux/amd64 \
  ${IMG_NAME}:${IMG_VERSION}
```

Push

```bash
REGION=us-east1
PROJECT=c-tc-429521
docker tag ${IMG_NAME}:${IMG_VERSION} \
  ${REGION}-docker.pkg.dev/${PROJECT}/${IMG_NAME}/${IMG_NAME}:${IMG_VERSION} \
  && docker push ${REGION}-docker.pkg.dev/${PROJECT}/${IMG_NAME}/${IMG_NAME}:${IMG_VERSION}
```

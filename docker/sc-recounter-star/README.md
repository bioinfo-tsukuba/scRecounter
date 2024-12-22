### sc-recounter-star

Build

> from the base directory of the repository

```bash
IMG_NAME=sc-recounter-star
IMG_VERSION=0.1.0
docker build \
  --file docker/${IMG_NAME}/Dockerfile \
  --build-arg CONDA_ENV_YAML=envs/star.yml \
  --platform linux/amd64 \
  --tag ${IMG_NAME}:${IMG_VERSION} \
  .
```

Push

```bash
REGION="us-east1"
PROJECT="c-tc-429521"
docker tag ${IMG_NAME}:${IMG_VERSION} \
  ${REGION}-docker.pkg.dev/${PROJECT}/${IMG_NAME}/${IMG_NAME}:${IMG_VERSION} \
  && docker push ${REGION}-docker.pkg.dev/${PROJECT}/${IMG_NAME}/${IMG_NAME}:${IMG_VERSION}
```

### auto-run

Env vars

```bash
IMG_NAME=auto-run
IMG_VERSION=0.1.0
source .env
```

Build the image

```bash
docker build \
  --file docker/${IMG_NAME}/Dockerfile \
  --build-arg GITHUB_PAT=${GITHUB_PAT} \
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
  --env DYNACONF="TEST" \
  --env GCP_SQL_DB_HOST=${GCP_SQL_DB_HOST} \
  --env GCP_SQL_DB_NAME=${GCP_SQL_DB_NAME} \
  --env GCP_SQL_DB_USERNAME=${GCP_SQL_DB_USERNAME} \
  --env GCP_SQL_DB_PASSWORD=${GCP_SQL_DB_PASSWORD} \
  --platform linux/amd64 \
  ${IMG_NAME}:${IMG_VERSION} \
  --max-records 10
```
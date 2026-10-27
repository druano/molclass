# Production (Erasmus MC)

On every new release a Github action will trigger a new Docker build for `ghcr.io/erasmusmc-bioinformatics/molclass`.  
The production environment uses this image in `docker-compose.yml` to run the Molclass production server.  
The production Molclass server is running at CIT, first ssh to the jumpoff server, `creig.erasmusmc.nl`, then ssh to `p-molclass01`.
Execute a `docker-compose pull` to pull the newly released image and then `docker-compose up -d --build` to restart the service with the updated release image.

## SSL certificate update

Make a copy of the old cert/key:  
```bash
cp /etc/nginx/certs/bifrost.erasmusmc.nl.crt /etc/nginx/certs/bifrost.erasmusmc.nl.crt.2023.old
cp /etc/nginx/certs/bifrost.erasmusmc.nl.key /etc/nginx/certs/bifrost.erasmusmc.nl.key.2023.old
```

Overwrite with the new cert/key:  
```bash
mv bifrost.erasmusmc.nl.crt /etc/nginx/certs/bifrost.erasmusmc.nl.crt
mv bifrost.erasmusmc.nl.key /etc/nginx/certs/bifrost.erasmusmc.nl.key
```

Restart the httpproxy docker:  
```bash
docker restart httpproxy
```

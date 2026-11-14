#!/usr/bin/env python3
"""Check Docker Hub's granted pull/push scope without uploading an image.

The token is read directly from Docker's HTTPS authorization endpoint and is
never logged or persisted. See https://docs.docker.com/reference/api/registry/auth/
"""
import base64
import json
import os
import urllib.parse
import urllib.request


def main():
    repository = "naotokubota/shiba"
    username = os.environ["DOCKER_HUB_USERNAME"]
    password = os.environ["DOCKER_HUB_ACCESS_TOKEN"]
    if not username or not password:
        raise RuntimeError("Docker Hub credentials are missing")
    query = urllib.parse.urlencode({"service": "registry.docker.io", "scope": f"repository:{repository}:pull,push"})
    authorization = base64.b64encode(f"{username}:{password}".encode()).decode()
    request = urllib.request.Request("https://auth.docker.io/token?" + query,
                                     headers={"Authorization": "Basic " + authorization})
    with urllib.request.urlopen(request, timeout=30) as response:
        token = json.load(response)["token"]
    # Inspect only the grant delivered over the authenticated HTTPS connection.
    # This is a preflight diagnostic, not a substitute for registry validation.
    encoded = token.split(".")[1]
    claims = json.loads(base64.urlsafe_b64decode(encoded + "=" * (-len(encoded) % 4)))
    grants = [entry for entry in claims.get("access", [])
              if entry.get("type") == "repository" and entry.get("name") == repository]
    if not any({"pull", "push"}.issubset(entry.get("actions", [])) for entry in grants):
        raise RuntimeError(f"Docker Hub did not grant pull and push access to {repository}")
    print(f"Docker Hub granted pull and push access to {repository}; no image was uploaded.")


if __name__ == "__main__":
    main()
